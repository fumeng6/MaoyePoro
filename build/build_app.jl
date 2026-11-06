using PackageCompiler

root = dirname(@__DIR__)
platform = Sys.iswindows() ? "Windows-x64" : Sys.islinux() ? "Linux-x86_64" : Sys.MACHINE
destination = isempty(ARGS) ? joinpath(root,"local_archive","builds","MaoyePoro-"*platform) : abspath(only(ARGS))

# PackageCompiler starts a separate Julia process for the workload; the build
# process's --threads option is not inherited by that process.
withenv("JULIA_NUM_THREADS" => string(max(2,Threads.nthreads()))) do
    create_app(root,destination;
        executables=["myporo"=>"julia_main","mycfs"=>"coulomb_main"],
        precompile_execution_file=joinpath(@__DIR__,"precompile_app.jl"),
        incremental=false,filter_stdlibs=false,cpu_target="generic",
        include_lazy_artifacts=false,include_transitive_dependencies=false)
end

if Sys.iswindows()
    powershell = joinpath(ENV["SystemRoot"],"System32","WindowsPowerShell","v1.0","powershell.exe")
    for executable in ("myporo.exe","mycfs.exe")
        run(`$powershell -NoProfile -ExecutionPolicy Bypass -File $(joinpath(@__DIR__,"add_manifest.ps1")) -Executable $(joinpath(destination,"bin",executable))`)
    end
end

for directory in ("cases","docs")
    cp(joinpath(root,directory),joinpath(destination,directory))
end
for file in ("README.md","LICENSE")
    cp(joinpath(root,file),joinpath(destination,file))
end
cp(normpath(Sys.BINDIR,"..","LICENSE.md"),joinpath(destination,"Julia-LICENSE.md"))
if Sys.iswindows()
    for file in ("myporo.cmd","mycfs.cmd","install.ps1","install.cmd")
        cp(joinpath(@__DIR__,file),joinpath(destination,file))
    end
end
if Sys.islinux()
    cp(joinpath(@__DIR__,"install.sh"),joinpath(destination,"install.sh"))
    chmod(joinpath(destination,"install.sh"),0o755)
end
println("Application bundle: ",destination)
