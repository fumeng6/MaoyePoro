# Coulomb failure stress, with tensors expressed in north/east/down coordinates.
# Stress is positive in tension; pore pressure is positive for an increase.

function fault_vectors(strike, dip, rake)
    st, di, ra = deg2rad.((strike, dip, rake))
    normal = [-sin(di)*sin(st), sin(di)*cos(st), -cos(di)]
    along_strike = [cos(st), sin(st), 0.0]
    down_dip = [-cos(di)*sin(st), cos(di)*cos(st), sin(di)]
    slip = cos(ra)*along_strike - sin(ra)*down_dip
    return normal, slip, along_strike, down_dip
end

"""
    coulomb_stress(stress, pressure, friction, strike, dip, rake)

Resolve a symmetric total stress tensor (Pa, tension positive, north/east/down)
onto a receiver fault. CFS = shear + friction*(normal_stress + pressure).
Rake is positive up dip; angles are in degrees. An undefined optimal rake is NaN.
"""
function coulomb_stress(stress, pressure, friction, strike, dip, rake)
    friction >= 0 || error("Friction must be nonnegative.")
    normal, slip, along_strike, down_dip = fault_vectors(strike, dip, rake)
    traction = stress * normal
    normal_stress = dot(normal, traction)
    shear = dot(slip, traction)
    shear_vector = traction - normal_stress*normal
    shear_max = norm(shear_vector)
    optimal_rake = shear_max <= eps(Float64)*norm(stress) ? NaN :
        mod(rad2deg(atan(-dot(down_dip, shear_vector), dot(along_strike, shear_vector))), 360)
    return (cfs = shear + friction*(normal_stress + pressure),
            normal_stress = normal_stress,
            shear = shear,
            cfs_optimal_rake = shear_max + friction*(normal_stress + pressure),
            optimal_rake = optimal_rake)
end

function mechanism_from_vectors(normal, slip)
    if normal[3] > 0
        normal, slip = -normal, -slip
    end
    strike = atan(normal[2], normal[1]) - pi/2
    dip = atan(hypot(normal[1], normal[2]), abs(normal[3]))
    along_strike = [cos(strike), sin(strike), 0.0]
    down_dip = [-cos(dip)*sin(strike), cos(dip)*cos(strike), sin(dip)]
    rake = atan(-dot(slip, down_dip), dot(slip, along_strike))
    return (strike = mod(rad2deg(strike), 360), dip = rad2deg(dip), rake = mod(rad2deg(rake), 360))
end

"""Maximum CFS over all planes and slip directions, from the principal stresses."""
function optimal_coulomb_stress(stress, pressure, friction; reference = (0.0, 90.0, 0.0))
    friction >= 0 || error("Friction must be nonnegative.")
    eig = eigen(Symmetric(stress))
    smin, smid, smax = eig.values
    delta = smax - smin
    normal_stress = 0.5*(smax+smin + delta*friction/hypot(1, friction))
    cfs = 0.5*delta*hypot(1, friction) + friction*(0.5*(smax+smin)+pressure)
    tol = 1e-10*norm(eig.values)
    if delta <= tol
        return (cfs = cfs, normal_stress = normal_stress, planes = nothing,
                orientations_unique = false)
    end
    alpha = 0.5*atan(1, friction)
    vmax, vmin = eig.vectors[:, 3], eig.vectors[:, 1]
    normals = (vmax*cos(alpha)-vmin*sin(alpha), vmax*cos(alpha)+vmin*sin(alpha))
    slips = (vmax*sin(alpha)+vmin*cos(alpha), vmax*sin(alpha)-vmin*cos(alpha))
    planes = mechanism_from_vectors.(normals, slips)
    nref, tref, _, _ = fault_vectors(reference...)
    correlations = [dot(nref, normals[i])*dot(tref, slips[i]) for i in 1:2]
    if correlations[2] > correlations[1]
        planes = reverse(planes)
    end
    return (cfs = cfs, normal_stress = normal_stress, planes = planes,
            orientations_unique = smax-smid > tol && smid-smin > tol)
end

"""Rotate a cylindrical (r,t,z) tensor to north/east/down; azimuth points away from the source."""
function cylindrical_to_ned(tensor, azimuth)
    c, s = cosd(azimuth), sind(azimuth)
    rotation = [c -s 0.0; s c 0.0; 0.0 0.0 1.0]
    return rotation * tensor * transpose(rotation)
end

function receiver_radial_azimuth(lat0, lon0, latr, lonr)
    source_lat, receiver_lat, dlon = deg2rad.((lat0, latr, lon0-lonr))
    east = sin(dlon)*cos(source_lat)
    north = cos(receiver_lat)*sin(source_lat) - sin(receiver_lat)*cos(source_lat)*cos(dlon)
    return mod(rad2deg(atan(east, north)) + 180, 360)
end

"""
    coulomb_from_strain(; lat0, lon0, latr, lonr, mu, nu, alpha, strain, pressure, friction, strike, dip, rake)

Form the coupled total stress λ tr(ε)I+2με-alpha*pressure*I,
then resolve it onto the fault and add +pressure in the Coulomb normal term.
`alpha` is the Biot coefficient at the receiver and must be supplied explicitly.
`strain` order is (ezz, err, ett, ezr, ert, etz); shear strains are tensor components.
Strain and pressure are changes from the same reference state. Both returned
stress tensors are in north/east/down coordinates (Pa, tension positive).
"""
function coulomb_from_strain(; lat0, lon0, latr, lonr, mu, nu, alpha, strain, pressure,
                  friction, strike, dip, rake)
    mu > 0 && 0 < nu < 0.5 || error("Require μ>0 and 0<ν<0.5.")
    0 <= alpha <= 1+8eps(Float64) || error("Require a receiver Biot coefficient 0≤alpha≤1.")
    ezz, err, ett, ezr, ert, etz = strain
    lambda = 2mu*nu/(1-2nu)
    epsilon = [err ert ezr; ert ett etz; ezr etz ezz]
    scyl = 2mu*epsilon + lambda*(err+ett+ezz)*I
    azimuth = receiver_radial_azimuth(lat0, lon0, latr, lonr)
    matrix_stress = cylindrical_to_ned(scyl, azimuth)
    total_stress = matrix_stress - alpha*pressure*I
    receiver = coulomb_stress(total_stress, pressure, friction, strike, dip, rake)
    optimal = optimal_coulomb_stress(total_stress, pressure, friction; reference = (strike, dip, rake))
    return (matrix_stress = matrix_stress, total_stress = total_stress,
            receiver = receiver, optimal = optimal,
            radial_azimuth = azimuth, convention = :total_stress_plus_pressure)
end

"""Read named Coulomb-stress parameters from TOML.

Coordinates and fault angles are in degrees; pressure and shear modulus are in Pa.
`strain` holds dimensionless cylindrical tensor components, including tensor shear
strains. Strain and pressure are changes from the same reference state.
"""
function read_coulomb_input(path::AbstractString)
    lowercase(splitext(path)[2]) == ".toml" ||
        throw(ArgumentError("Coulomb stress input must be a .toml file."))
    config = TOML.parsefile(path)
    check_config_keys(config, ("schema_version", "pressure_pa", "source", "receiver",
        "material", "fault", "strain"), "Coulomb root")
    version = get(config, "schema_version", 1)
    version isa Integer && version == 1 ||
        throw(ArgumentError("Unsupported Coulomb input schema_version: $version."))
    source = check_config_keys(config["source"], ("latitude_deg", "longitude_deg"), "source")
    receiver = check_config_keys(config["receiver"], ("latitude_deg", "longitude_deg"), "receiver")
    material = check_config_keys(config["material"],
        ("shear_modulus_pa", "poisson", "biot_alpha"), "material")
    fault = check_config_keys(config["fault"],
        ("strike_deg", "dip_deg", "rake_deg", "friction"), "fault")
    strain = check_config_keys(config["strain"], ("ezz", "err", "ett", "ezr", "ert", "etz"), "strain")
    return (lat0=Float64(source["latitude_deg"]), lon0=Float64(source["longitude_deg"]),
            latr=Float64(receiver["latitude_deg"]), lonr=Float64(receiver["longitude_deg"]),
            mu=Float64(material["shear_modulus_pa"]), nu=Float64(material["poisson"]),
            alpha=Float64(material["biot_alpha"]), pressure=Float64(config["pressure_pa"]),
            strike=Float64(fault["strike_deg"]), dip=Float64(fault["dip_deg"]),
            rake=Float64(fault["rake_deg"]), friction=Float64(fault["friction"]),
            strain=Tuple(Float64(strain[name]) for name in ("ezz", "err", "ett", "ezr", "ert", "etz")))
end

"""Compute Coulomb stress from a strain input file and print the numerical result."""
function run_coulomb(path::AbstractString; io = stdout)
    result = coulomb_from_strain(; read_coulomb_input(path)...)
    println(io, "MaoyePoro Coulomb stress: coupled total stress from strain (including -alpha*p*I), CFS = shear + friction*(normal stress + pore pressure) (Pa; tension positive)")
    s = result.total_stress
    @printf(io, "Snn See Sdd Sne Sed Sdn [Pa]\n%.10e %.10e %.10e %.10e %.10e %.10e\n",
            s[1,1], s[2,2], s[3,3], s[1,2], s[2,3], s[3,1])
    r = result.receiver
    @printf(io, "CFS_Rec SIG_Rec CFS_Rec_Opt [Pa] RAKE_Rec_Opt [deg]\n%.10e %.10e %.10e %.8f\n",
            r.cfs, r.normal_stress, r.cfs_optimal_rake, r.optimal_rake)
    @printf(io, "CFS_3D_Opt SIG_3D_Opt [Pa]\n%.10e %.10e\n", result.optimal.cfs, result.optimal.normal_stress)
    if result.optimal.planes !== nothing
        println(io, "Optimal strike dip rake [deg] (nonunique if a principal stress repeats)")
        for plane in result.optimal.planes
            @printf(io, "%.8f %.8f %.8f\n", plane.strike, plane.dip, plane.rake)
        end
    else
        println(io, "Isotropic stress: optimal plane orientation is undefined.")
    end
    return result
end
