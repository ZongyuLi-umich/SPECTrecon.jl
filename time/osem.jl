# project.jl

using BenchmarkTools: @btime
using SPECTrecon
using LinearMapsAA
using SPECTrecon: osem!
using MATLAB

function call_SPECTosem_matlab(mpath, proj, scatter, mumap, psfs, dy, nthreads)

    mat"""
    addpath($mpath)
    SPECTosem_matlab($proj, $scatter, $mumap, $psfs, $dy, $nthreads);
    """
end

function osem_time()
    T = Float32
    nx = 128
    ny = 128
    nz = 80
    nview = 128

    mumap = rand(T, nx, ny, nz)

    nx_psf = 37
    nz_psf = 37
    psfs = rand(T, nx_psf, nz_psf, ny, nview)
    psfs = psfs .+ mapslices(reverse, psfs, dims = [1, 2])
    psfs = psfs .+ mapslices(transpose, psfs, dims = [1, 2])
    psfs = psfs ./ mapslices(sum, psfs, dims = [1, 2])


    dy = T(4.7952)
    niter = 16
    nblocks = 4
    proj = rand(T, nx, nz, nview)
    scatter = rand(T, nx, nz, nview)

    plan = SPECTplan(mumap, psfs, dy)
    forw! = (y,x) -> project!(y, x, plan)
    back! = (x,y) -> backproject!(x, y, plan)
    idim = (nx,ny,nz)
    odim = (nx,nz,nview)
    A = LinearMapAA(forw!, back!, (prod(odim),prod(idim)); T, odim, idim)
    Ab = Ablock(plan, nblocks) # create a linear map for each block
    x0 = rand(T, nx, ny, nz)
    xhat = similar(x0)
    @time osem!(xhat, x0, proj, scatter, Ab; niter)

    mpath = pwd()
    println("osem-matlab")
    println("Warning: Check if MIRT is installed")
    nthreads = Threads.nthreads()
    call_SPECTosem_matlab(mpath, proj, scatter, mumap, psfs, dy, nthreads) # 216.518 ms, about 0.01 GiB
    nothing
end

# run all functions, time may vary on different machines, will allocate ~4MB memory.
osem_time()
osem_time()
#= one, fast
    378.348 ms (25983 allocations: 1.37 MiB)
one, mem
    793.918 ms (31257 allocations: 1.76 MiB)
two, fast
    289.804 ms (25982 allocations: 1.37 MiB)
two, mem
    618.369 ms (31268 allocations: 1.76 MiB)
MIRT
    230.125 ms
=#
