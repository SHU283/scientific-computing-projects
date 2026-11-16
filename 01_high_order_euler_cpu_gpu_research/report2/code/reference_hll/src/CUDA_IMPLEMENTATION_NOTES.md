# CUDA execution path

The command

```text
./final_project --backend=cuda --method=weno2 --problem=gaussian_advection_2d
```

follows this high-level path:

```text
main.cpp
  -> run_single_case(config)
  -> select the CUDA backend
  -> cuda_backend_probe()
  -> construct the grid and initial state
  -> cuda_run_weno2_2d(config, grid, initial_state)

cuda_backend.cu
  -> allocate device arrays with cudaMalloc
  -> copy the initial state from host to device
  -> for each time step:
       fp_cuda_weno2_rhs_2d_kernel<<<...>>>
       fp_cuda_add_scaled_2d_kernel<<<...>>>
       fp_cuda_weno2_rhs_2d_kernel<<<...>>>
       fp_cuda_rk2_finalize_2d_kernel<<<...>>>
       cudaDeviceSynchronize()
  -> copy the final state from device to host
  -> return SimulationResult2D
```

The evolving field remains resident in device memory throughout the time loop. Host-device transfers therefore occur at setup and final collection rather than after every numerical stage.
