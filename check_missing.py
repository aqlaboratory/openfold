#!/usr/bin/env python3

def check_missing_packages():
    """Check which packages from environment.yml are missing"""
    print("=" * 60)
    print("Checking packages from environment.yml")
    print("=" * 60)
    
    # Python packages that should be available
    python_packages = [
        'biopython', 'numpy', 'pandas', 'requests', 'scipy', 
        'tqdm', 'wandb', 'openmm', 'pdbfixer', 'torch',
        'pytorch_lightning', 'deepspeed', 'flash_attn'
    ]
    
    print("Python packages:")
    print("-" * 30)
    missing_python = []
    for pkg in python_packages:
        try:
            __import__(pkg)
            print(f"✅ {pkg}")
        except ImportError as e:
            print(f"❌ {pkg} - NOT INSTALLED")
            missing_python.append(pkg)
    
    # Check system tools
    import subprocess
    import os
    
    print("\nSystem tools:")
    print("-" * 30)
    system_tools = ['aria2c', 'aws', 'git', 'hmmsearch', 'hhblits', 'kalign']
    missing_tools = []
    
    for tool in system_tools:
        try:
            result = subprocess.run(['which', tool], 
                                   capture_output=True, text=True)
            if result.returncode == 0:
                print(f"✅ {tool}")
            else:
                print(f"❌ {tool} - NOT INSTALLED")
                missing_tools.append(tool)
        except:
            print(f"❌ {tool} - NOT INSTALLED")
            missing_tools.append(tool)
    
    # Check CUDA tools
    print("\nCUDA tools:")
    print("-" * 30)
    cuda_tools = ['nvcc', 'nvidia-smi']
    missing_cuda = []
    
    for tool in cuda_tools:
        try:
            result = subprocess.run(['which', tool], 
                                   capture_output=True, text=True)
            if result.returncode == 0:
                print(f"✅ {tool}")
            else:
                print(f"❌ {tool} - NOT INSTALLED")
                missing_cuda.append(tool)
        except:
            print(f"❌ {tool} - NOT INSTALLED")
            missing_cuda.append(tool)
    
    print("\n" + "=" * 60)
    print("Summary of missing components:")
    print("=" * 60)
    
    if missing_python:
        print(f"Missing Python packages: {', '.join(missing_python)}")
    
    if missing_tools:
        print(f"Missing system tools: {', '.join(missing_tools)}")
        
    if missing_cuda:
        print(f"Missing CUDA tools: {', '.join(missing_cuda)}")
    
    if not missing_python and not missing_tools and not missing_cuda:
        print("🎉 All packages and tools are installed!")
    else:
        print("\nRecommendations:")
        if missing_tools:
            print("- Install bioinformatics tools with conda/mamba or apt")
        if missing_cuda:
            print("- CUDA tools missing but PyTorch CUDA works fine")

if __name__ == "__main__":
    check_missing_packages() 