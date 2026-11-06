using MaoyePoro, LinearAlgebra
const MP = MaoyePoro
const ROOT = dirname(@__DIR__)
BLAS.set_num_threads(1)
println("Precompile workload: ",Threads.nthreads()," Julia threads")
Threads.nthreads() >= 2 || error("Run the precompile workload with at least two Julia threads to compile parallel frequency evaluation.")

# Match the C executable launcher's parser entry, rather than Meta.parse(String).
ccall(:jl_eval_string,Any,(Cstring,),"nothing")
precompile(Tuple{typeof(Base.setproperty!),Module,Symbol,String})

MP.cli_main(["--help"])
MP.cli_main(["--version"])
MP.coulomb_cli(["--help"])
MP.coulomb_cli(["--version"])

mktempdir() do directory
    cd(directory) do
        MP.cli_main(["--template"])
        MP.cli_main(["--template","custom.toml"])
    end

    # Small nonzero histories compile time inversion, serial/threaded frequency
    # evaluation, all output formats, snapshots and the source-pressure integral.
    c = read_input(joinpath(ROOT,"cases","constant_injection.toml"))
    c = PoroelasticInput(c;r=[300.0,600.0],zr=[1400.0,1600.0],nt=9,twindow=1e6,
        ts0=[0.0,100000.0,1000000.0],sinj0=[0.0,0.01,0.01],
        timesn=[500000.0,562500.0],filesn=["saved.dat","interpolated.dat"],
        compute_equivalent=true)
    input_path = joinpath(directory,"case.toml")
    write_input(input_path,c;history_file="history.csv")
    write_input(joinpath(directory,"inline.toml"),c)
    MP.cli_main([input_path,"--check","--max-frequencies","512","--image-factor","0.001"])
    output_dir = joinpath(directory,"output")
    MP.cli_main([input_path,"--output",output_dir,"--format","both",
        "--max-frequencies","512","--quiet"])
    archive = joinpath(output_dir,"result.jld2")
    result = load_result(archive)
    read_field(archive,:pp)
    read_field(archive,"uz")
    read_timeseries(joinpath(output_dir,"pp.dat"))
    field(result,:pp)
    field(result,"uz")
    snapshot(result,500000.0)
    snapshot(result,562500.0)
    sample_field(result,:pp;r=300.0,z=1400.0,times=last(result.time))
    sample_field(archive,:pp;r=[300.0,600.0],z=[1400.0,1600.0],tspan=(500000.0,1000000.0))
    sample_field(result,"uz";r=[450.0],z=[1500.0],times=[562500.0],method=:linear)
    sample_field(archive,:pp;x=[180.0 360.0],y=[240.0 480.0],z=1400.0,times=500000.0)
    strain_tensor(result;x=180.0,y=240.0,z=1400.0,times=500000.0)
    strain_tensor(archive;r=[300.0,600.0],z=1400.0,azimuth_deg=[0.0,90.0],times=[500000.0])
    write_outputs(result;output_dir=joinpath(directory,"text"),format=:text)
    write_outputs(result;output_dir=joinpath(directory,"jld2"),format=:jld2)
    save_result(joinpath(directory,"compressed.jld2"),result;compress=true)
    load_result(joinpath(directory,"compressed.jld2"))
    MP.cli_main([input_path,"--output",joinpath(directory,"serial"),"--format","jld2",
        "--max-frequencies","512","--serial","--quiet"])
    MP.cli_main([input_path])
    MP.cli_main([input_path,"--output",joinpath(directory,"default")])
    MP.cli_main([input_path,"--output",joinpath(directory,"default"),"--quiet"])
    MP.cli_main([input_path,"--output",joinpath(directory,"default"),"--format","both"])
    MP.cli_main([input_path,"--output",joinpath(directory,"default"),"--format","both","--quiet"])
    MP.cli_main([input_path,"--output",joinpath(directory,"serial_text"),"--format","text","--serial","--quiet"])

    for name in ("volume_change","pressure_moment")
        println("Precompiling nonzero ",name)
        source = PoroelasticInput(c;istype=name == "volume_change" ? 1 : 3,
            zstop=1000.0,zsbtm=1100.0,r=[300.0],zr=[1500.0],nt=5,
            ts0=[0.0,1e6],sinj0=name == "volume_change" ? [0.0,1000.0] : [0.0,1e10],
            timesn=Float64[],filesn=String[])
        run_simulation(source;output_dir=joinpath(directory,name),format=:both,
            verbose=false,max_frequencies=512,parallel=false)
        load_result(joinpath(directory,name,"result.jld2"))
    end

    println("Precompiling nonzero surface force and drained elastic limit")
    force = PoroelasticInput(c;istype=0,zstop=0.0,zsbtm=0.0,sradius=500.0,
        isurfcon=1,r=[0.0,1000.0],zr=[0.0],nt=5,
        ts0=[0.0,1e6],sinj0=[0.0,1e9],
        timesn=Float64[],filesn=String[],compute_equivalent=false)
    model = build_model(force)
    ws = MP.KernelWorkspace(model,build_grid(force,model),force)
    MP.field_spectrum(ws,1e-6+2e-6im,force.r,[0.0],0,
        [MP.bessel_table(force.r,force.sradius,0.001;nbsj=256)];extend=false)
    elastic_nodes = copy(force.model_nodes)
    elastic_nodes[:,4] .= elastic_nodes[:,3]
    elastic = PoroelasticInput(force;model_nodes=elastic_nodes)
    run_simulation(elastic;output_dir=joinpath(directory,"elastic"),
        format=:both,verbose=false,max_frequencies=512)

    # Compile the real layered/half-space kernels at nonzero complex frequencies.
    # These are precompile workloads, not numerical-convergence comparisons.
    # Avoid repeatedly inverting a full layered time history during every build.
    layered = read_input(joinpath(ROOT,"cases","layered_injection.toml"))
    for surface in (0,1,2)
        println("Precompiling layered kernels, surface condition ",surface)
        source = PoroelasticInput(layered;isurfcon=surface,r=[400.0],zr=[0.0,1400.0],nt=5,
            timesn=Float64[],filesn=String[])
        model = build_model(source)
        grid = build_grid(source,model)
        ws = MP.KernelWorkspace(model,grid,source)
        zrs2 = [min((z-source.zstop)^2,(z-source.zsbtm)^2) for z in source.zr]
        table = MP.bessel_table(source.r,source.sradius,0.001;nbsj=256)
        for s in (1e-6+2e-6im,0.01+0.02im)
            sf = MP.source_vector(ws,s)
            pilot,_ = MP.field_spectrum(ws,s,source.r,zrs2,0,[table];extend=false)
            atol = max.(abs.(pilot),1e-12)
            MP.field_spectrum(ws,s,source.r,zrs2,0,[table];atol,extend=false)
            MP.volume_direct_spectrum(ws,s,source.r,sf,atol)
            MP.mapped_gauss_correction(ws,s,sf,source.r,atol)
        end
    end

    # Near-source integration uses the Hankel path instead of separated-volume
    # quadrature; cover the source-normalization branch for prescribed pressure.
    near = PoroelasticInput(c;istype=3,zstop=1000.0,zsbtm=1100.0,
        r=[50.0],zr=[1110.0],ts0=[0.0,1e6],sinj0=[0.0,1e10],
        timesn=Float64[],filesn=String[])
    model = build_model(near)
    ws = MP.KernelWorkspace(model,build_grid(near,model),near)
    MP.field_spectrum(ws,0.01+0.02im,near.r,[100.0],0,
        [MP.bessel_table(near.r,near.sradius,0.001;nbsj=256)];extend=false)

    # C2: exterior receivers closer than one source radius use the fluid-only
    # volume quadrature, rather than the separated all-field volume path.
    exterior = PoroelasticInput(near;zstop=900.0,zsbtm=1100.0,sradius=100.0,
        r=[150.0],zr=[900.0,1000.0,1100.0],isurfcon=0,
        twindow=2000.0,nt=9,ts0=[0.0,1000.0,2000.0],sinj0=[0.0,1e10,1e10],
        model_nodes=[0.0 1.5e10 0.25 0.36 0.5 0.024],
        timesn=Float64[],filesn=String[],compute_equivalent=false)
    run_simulation(exterior;output_dir=joinpath(directory,"exterior_pressure"),
        format=:both,verbose=false,max_frequencies=512)

    # T3: a buried force disk, a sealed surface, a material interface and a
    # true history jump. Compile spectra and the t=0 half derivative directly;
    # a build workload need not repeat the full T3 convergence calculation.
    buried = PoroelasticInput(exterior;istype=0,zstop=1300.0,zsbtm=1300.0,
        sradius=200.0,r=[0.0,300.0],zr=[0.0,700.0,1300.0,2000.0],
        isurfcon=2,twindow=1e5,nt=5,ts0=[0.0,0.0,1e5],sinj0=[0.0,1e9,1e9],
        model_nodes=[0.0 2e10 .25 .35 .5 .02;700.0 2e10 .25 .35 .5 .02;
                     700.0 5e9 .25 .36 .5 1.0])
    model=build_model(buried)
    ws=MP.KernelWorkspace(model,build_grid(buried,model),buried)
    zrs2=(buried.zr .- buried.zstop).^2
    tables=[MP.bessel_table(buried.r,buried.sradius,0.001;nbsj=512)]
    df=0.5/buried.twindow
    abscissa=-df*log(0.001)
    spectra=[first(MP.field_spectrum(ws,complex(abscissa,2pi*j*df),
        buried.r,zrs2,0,tables;extend=false)) for j in 0:1]
    _,coeff,_=MP.field_asymptotics(ws,buried.r,zrs2,0,tables,abs.(spectra[1]))
    MP.inverse_fields(spectra,df,abscissa,buried,coeff)

    stress_path=joinpath(directory,"stress.toml")
    write(stress_path,raw"""
    schema_version = 1
    pressure_pa = 1000.0
    [source]
    latitude_deg = 30.0
    longitude_deg = 103.0
    [receiver]
    latitude_deg = 30.002
    longitude_deg = 103.003
    [material]
    shear_modulus_pa = 1.0e10
    poisson = 0.25
    biot_alpha = 0.8
    [fault]
    strike_deg = 30.0
    dip_deg = 60.0
    rake_deg = 90.0
    friction = 0.6
    [strain]
    ezz = 1.0e-6
    err = -2.0e-6
    ett = 0.5e-6
    ezr = 0.3e-6
    ert = 0.0
    etz = 0.0
    """)
    MP.coulomb_cli([stress_path])
    open(joinpath(directory,"coulomb.log"),"w") do io
        # Keep a standalone specialization instead of inlining it into this workload.
        Base.invokelatest(run_coulomb,stress_path;io)
    end
    Base.invokelatest(run_coulomb,stress_path;io=IOContext(stdout))
    MP.optimal_coulomb_stress(zeros(3,3),0.0,0.6)
end
