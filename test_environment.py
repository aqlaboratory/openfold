#!/usr/bin/env python3
import sys

def test_pytorch():
    """测试PyTorch和CUDA"""
    try:
        import torch
        print(f'✅ PyTorch version: {torch.__version__}')
        print(f'✅ CUDA available: {torch.cuda.is_available()}')
        if torch.cuda.is_available():
            print(f'✅ CUDA device: {torch.cuda.get_device_name()}')
            print(f'✅ CUDA version: {torch.version.cuda}')
            print(f'✅ CUDA device count: {torch.cuda.device_count()}')
        return True
    except Exception as e:
        print(f'❌ PyTorch test failed: {e}')
        return False

def test_openfold():
    """测试OpenFold导入"""
    try:
        import openfold
        print('✅ OpenFold package imported successfully')
        
        # 测试核心模块导入
        from openfold.model.model import AlphaFold
        print('✅ AlphaFold model imported successfully')
        
        from openfold.data import data_pipeline
        print('✅ Data pipeline imported successfully')
        
        from openfold.utils.kernel import attention_core
        print('✅ CUDA attention kernel imported successfully')
        
        return True
    except Exception as e:
        print(f'❌ OpenFold test failed: {e}')
        return False

def test_dependencies():
    """测试关键依赖"""
    dependencies = [
        'numpy', 'pandas', 'scipy', 'ml_collections',
        'pytorch_lightning', 'biopython', 'openmm', 'pdbfixer'
    ]
    
    success = True
    for dep in dependencies:
        try:
            if dep == 'pytorch_lightning':
                import lightning
                print(f'✅ {dep} (lightning) imported')
            elif dep == 'biopython':
                import Bio
                print(f'✅ {dep} imported')
            elif dep == 'ml_collections':
                import ml_collections
                print(f'✅ {dep} imported')
            else:
                __import__(dep)
                print(f'✅ {dep} imported')
        except ImportError as e:
            print(f'❌ {dep} import failed: {e}')
            success = False
    
    return success

def test_cuda_extensions():
    """测试CUDA扩展"""
    try:
        import attn_core_inplace_cuda
        print('✅ CUDA attention extension loaded successfully')
        return True
    except Exception as e:
        print(f'❌ CUDA extension test failed: {e}')
        return False

if __name__ == "__main__":
    print("=== OpenFold环境测试 ===\n")
    
    print("1. 测试PyTorch和CUDA:")
    pytorch_ok = test_pytorch()
    print()
    
    print("2. 测试依赖库:")
    deps_ok = test_dependencies()
    print()
    
    print("3. 测试CUDA扩展:")
    cuda_ok = test_cuda_extensions()
    print()
    
    print("4. 测试OpenFold:")
    openfold_ok = test_openfold()
    print()
    
    print("=== 测试总结 ===")
    if all([pytorch_ok, deps_ok, cuda_ok, openfold_ok]):
        print("🎉 所有测试通过！OpenFold环境配置成功")
        print("📝 注意：仍需下载模型参数文件才能进行结构预测")
    else:
        print("⚠️ 部分测试失败，请检查上述错误信息")
        if not pytorch_ok:
            print("   - PyTorch/CUDA 配置有问题")
        if not deps_ok:
            print("   - 依赖库缺失")
        if not cuda_ok:
            print("   - CUDA扩展有问题")
        if not openfold_ok:
            print("   - OpenFold安装有问题") 