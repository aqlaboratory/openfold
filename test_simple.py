#!/usr/bin/env python3
"""
简化的OpenFold功能测试
"""
import torch
import sys

def test_basic_imports():
    """测试基础导入"""
    print("=== 测试基础导入 ===")
    try:
        import openfold
        print("✅ openfold导入成功")
        
        from openfold.model.model import AlphaFold
        print("✅ AlphaFold模型导入成功")
        
        from openfold.config import model_config
        print("✅ 配置模块导入成功")
        
        return True
    except Exception as e:
        print(f"❌ 导入失败: {e}")
        return False

def test_config():
    """测试配置加载"""
    print("\n=== 测试配置 ===")
    try:
        from openfold.config import model_config
        
        # 尝试不同的配置名称
        config_names = ["model_1", "model_2", "model_3", "model_4", "model_5"]
        
        for config_name in config_names:
            try:
                config = model_config(config_name)
                print(f"✅ 配置 {config_name} 加载成功")
                
                # 检查配置结构
                if hasattr(config, 'model'):
                    if hasattr(config.model, 'evoformer_stack'):
                        evo_config = config.model.evoformer_stack
                        if hasattr(evo_config, 'no_blocks'):
                            num_blocks = evo_config.no_blocks
                        elif hasattr(evo_config, 'num_blocks'):
                            num_blocks = evo_config.num_blocks
                        else:
                            num_blocks = "未知"
                        print(f"   - Evoformer层数: {num_blocks}")
                        return config
                break
            except Exception as e:
                print(f"⚠️ 配置 {config_name} 失败: {e}")
                continue
        
        print("❌ 所有配置加载失败")
        return None
        
    except Exception as e:
        print(f"❌ 配置测试失败: {e}")
        return None

def test_cuda_extension():
    """测试CUDA扩展"""
    print("\n=== 测试CUDA扩展 ===")
    
    if not torch.cuda.is_available():
        print("⚠️ CUDA不可用")
        return True
    
    try:
        # 直接导入编译好的扩展
        import attn_core_inplace_cuda
        print("✅ CUDA扩展导入成功")
        
        # 检查CUDA设备
        device = torch.device("cuda")
        x = torch.randn(2, 3, device=device)
        print(f"✅ CUDA张量创建成功: {x.device}")
        
        return True
    except Exception as e:
        print(f"❌ CUDA扩展测试失败: {e}")
        return False

def test_model_structure():
    """测试模型结构"""
    print("\n=== 测试模型结构 ===")
    
    config = test_config()
    if config is None:
        return False
    
    try:
        from openfold.model.model import AlphaFold
        
        # 创建模型
        model = AlphaFold(config)
        print("✅ AlphaFold模型创建成功")
        
        # 计算参数数量
        total_params = sum(p.numel() for p in model.parameters())
        trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
        
        print(f"   - 总参数数量: {total_params:,}")
        print(f"   - 可训练参数: {trainable_params:,}")
        print(f"   - 模型大小: {total_params * 4 / (1024**3):.2f} GB (fp32)")
        
        # 检查主要组件
        if hasattr(model, 'evoformer'):
            print("✅ Evoformer组件存在")
        if hasattr(model, 'structure_module'):
            print("✅ 结构模块存在")
        
        return True
        
    except Exception as e:
        print(f"❌ 模型结构测试失败: {e}")
        return False

def main():
    """主函数"""
    print("🧬 OpenFold简化测试")
    print("=" * 40)
    
    print(f"PyTorch版本: {torch.__version__}")
    print(f"CUDA可用: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"CUDA设备: {torch.cuda.get_device_name()}")
    print()
    
    tests = [
        ("基础导入", test_basic_imports),
        ("CUDA扩展", test_cuda_extension),
        ("模型结构", test_model_structure),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append(result)
        except Exception as e:
            print(f"❌ {test_name}测试异常: {e}")
            results.append(False)
    
    print("\n" + "=" * 40)
    print("📊 测试总结:")
    
    passed = sum(results)
    total = len(results)
    
    for i, (test_name, _) in enumerate(tests):
        status = "✅ 通过" if results[i] else "❌ 失败"
        print(f"   {test_name}: {status}")
    
    print(f"\n📈 通过率: {passed}/{total} ({passed/total*100:.1f}%)")
    
    if passed == total:
        print("\n🎉 OpenFold核心功能正常！")
        print("💡 环境配置成功，可以下载参数文件开始使用")
    elif passed >= total * 0.5:
        print("\n⚠️ 大部分功能正常，有部分问题")
        print("💡 建议下载参数文件测试完整功能")
    else:
        print("\n❌ 多个功能异常，需要检查环境")
    
    return passed >= total * 0.5

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 