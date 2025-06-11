#!/usr/bin/env python3
"""
测试OpenFold参数加载
"""
import torch
import os
import sys

def test_alphafold_params():
    """测试AlphaFold参数"""
    print("=== 测试AlphaFold2参数 ===")
    
    params_dir = "openfold/resources/params"
    if not os.path.exists(params_dir):
        print("❌ AlphaFold参数目录不存在")
        return False
    
    # 检查主要模型文件
    model_files = [
        "params_model_1.npz",
        "params_model_2.npz", 
        "params_model_3.npz",
        "params_model_4.npz",
        "params_model_5.npz"
    ]
    
    found_models = []
    for model_file in model_files:
        full_path = os.path.join(params_dir, model_file)
        if os.path.exists(full_path):
            size_mb = os.path.getsize(full_path) / (1024*1024)
            print(f"✅ {model_file}: {size_mb:.1f}MB")
            found_models.append(model_file)
        else:
            print(f"❌ {model_file}: 缺失")
    
    print(f"AlphaFold模型: {len(found_models)}/5 个可用")
    return len(found_models) >= 3  # 至少需要3个模型

def test_openfold_params():
    """测试OpenFold参数"""
    print("\n=== 测试OpenFold参数 ===")
    
    params_dir = "openfold/resources/openfold_params"
    if not os.path.exists(params_dir):
        print("❌ OpenFold参数目录不存在")
        return False
    
    # 检查所有.pt文件
    pt_files = [f for f in os.listdir(params_dir) if f.endswith('.pt')]
    pt_files.sort()
    
    if not pt_files:
        print("❌ 没有找到.pt参数文件")
        return False
    
    print(f"找到 {len(pt_files)} 个参数文件:")
    total_size = 0
    valid_files = 0
    
    for pt_file in pt_files:
        full_path = os.path.join(params_dir, pt_file)
        size_mb = os.path.getsize(full_path) / (1024*1024)
        total_size += size_mb
        
        # 尝试加载文件头部验证
        try:
            with open(full_path, 'rb') as f:
                # 检查是否是有效的PyTorch文件
                header = f.read(8)
                if header[:2] == b'PK':  # ZIP格式（PyTorch使用）
                    status = "✅"
                    valid_files += 1
                else:
                    status = "⚠️"
            print(f"   {status} {pt_file}: {size_mb:.1f}MB")
        except Exception as e:
            print(f"   ❌ {pt_file}: 无法读取 ({e})")
    
    print(f"总大小: {total_size:.1f}MB ({total_size/1024:.2f}GB)")
    print(f"有效文件: {valid_files}/{len(pt_files)}")
    
    return valid_files >= 3  # 至少需要3个有效模型

def test_model_loading():
    """测试模型加载"""
    print("\n=== 测试模型加载 ===")
    
    try:
        from openfold.config import model_config
        from openfold.model.model import AlphaFold
        
        # 创建模型
        config = model_config("model_1")
        model = AlphaFold(config)
        print("✅ 模型结构创建成功")
        
        # 检查参数文件是否可以被模型识别
        params_dir = "openfold/resources/openfold_params"
        pt_files = [f for f in os.listdir(params_dir) if f.endswith('.pt')]
        
        if pt_files:
            # 尝试加载一个参数文件（不实际加载到模型，只验证格式）
            test_file = os.path.join(params_dir, pt_files[0])
            try:
                checkpoint = torch.load(test_file, map_location='cpu')
                if isinstance(checkpoint, dict):
                    print(f"✅ 参数文件格式正确: {pt_files[0]}")
                    if 'ema' in checkpoint:
                        print("   - 包含EMA权重")
                    if 'model' in checkpoint:
                        print("   - 包含模型权重")
                    return True
                else:
                    print(f"⚠️ 参数文件格式异常: {pt_files[0]}")
                    return False
            except Exception as e:
                print(f"❌ 参数文件加载失败: {e}")
                return False
        else:
            print("❌ 没有找到参数文件")
            return False
            
    except Exception as e:
        print(f"❌ 模型加载测试失败: {e}")
        return False

def main():
    """主测试函数"""
    print("🧬 OpenFold参数完整性测试")
    print("=" * 40)
    
    tests = [
        ("AlphaFold2参数", test_alphafold_params),
        ("OpenFold参数", test_openfold_params),
        ("模型加载", test_model_loading),
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
    
    if passed == total:
        print(f"\n🎉 所有测试通过！OpenFold完全可用")
        print("💡 你可以开始进行蛋白质结构预测了！")
    elif passed >= 2:
        print(f"\n✅ 主要功能正常 ({passed}/{total})")
        print("💡 OpenFold基本可用")
    else:
        print(f"\n⚠️ 多个问题需要解决 ({passed}/{total})")
    
    return passed >= 2

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 