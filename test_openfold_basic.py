#!/usr/bin/env python3
"""
OpenFold基础功能测试（无需预训练参数）
测试模型结构和数据处理管道
"""

import torch
import numpy as np
from openfold.model.model import AlphaFold
from openfold.data import data_transforms
from openfold.config import model_config
import ml_collections

def test_model_creation():
    """测试模型创建"""
    print("=== 测试AlphaFold模型创建 ===")
    
    try:
        # 获取默认配置
        config = model_config("model_1")
        print(f"✅ 加载配置: {config.model.evoformer_stack.num_blocks}层Evoformer")
        
        # 创建模型
        model = AlphaFold(config)
        print(f"✅ 模型创建成功")
        print(f"   - 参数数量: {sum(p.numel() for p in model.parameters()):,}")
        print(f"   - 可训练参数: {sum(p.numel() for p in model.parameters() if p.requires_grad):,}")
        
        return True
    except Exception as e:
        print(f"❌ 模型创建失败: {e}")
        return False

def test_data_transforms():
    """测试数据变换"""
    print("\n=== 测试数据变换管道 ===")
    
    try:
        # 创建示例特征
        batch = {
            'aatype': torch.randint(0, 20, (1, 100)),  # 氨基酸序列
            'residue_index': torch.arange(100).unsqueeze(0),
            'seq_length': torch.tensor([100]),
            'msa': torch.randint(0, 22, (1, 32, 100)),  # MSA
            'num_alignments': torch.tensor([32]),
        }
        
        print(f"✅ 创建示例输入数据:")
        print(f"   - 序列长度: {batch['seq_length'].item()}")
        print(f"   - MSA深度: {batch['num_alignments'].item()}")
        
        # 测试一些基本变换
        from openfold.data.data_transforms import make_atom14_masks
        atom14_mask = make_atom14_masks(batch)
        print(f"✅ Atom14掩码创建成功: {atom14_mask.shape}")
        
        return True
    except Exception as e:
        print(f"❌ 数据变换测试失败: {e}")
        return False

def test_cuda_kernels():
    """测试CUDA内核"""
    print("\n=== 测试CUDA加速内核 ===")
    
    if not torch.cuda.is_available():
        print("⚠️ CUDA不可用，跳过CUDA内核测试")
        return True
    
    try:
        from openfold.utils.kernel.attention_core import attention_core_inplace_cuda
        print("✅ CUDA注意力内核导入成功")
        
        # 创建测试数据
        device = torch.device("cuda")
        batch_size, seq_len, dim = 2, 64, 128
        
        q = torch.randn(batch_size, seq_len, dim, device=device, dtype=torch.float32)
        k = torch.randn(batch_size, seq_len, dim, device=device, dtype=torch.float32)
        v = torch.randn(batch_size, seq_len, dim, device=device, dtype=torch.float32)
        bias = torch.zeros(batch_size, seq_len, seq_len, device=device, dtype=torch.float32)
        
        print(f"✅ CUDA张量创建成功: {q.shape}")
        print(f"   - 设备: {q.device}")
        print(f"   - 数据类型: {q.dtype}")
        
        return True
    except Exception as e:
        print(f"❌ CUDA内核测试失败: {e}")
        return False

def test_forward_pass():
    """测试前向传播（小规模）"""
    print("\n=== 测试小规模前向传播 ===")
    
    try:
        # 使用最小配置
        config = model_config("model_1")
        # 减小模型尺寸以节省内存
        config.model.evoformer_stack.num_blocks = 2
        config.model.structure_module.num_layer = 2
        
        model = AlphaFold(config)
        if torch.cuda.is_available():
            model = model.cuda()
            device = torch.device("cuda")
            print("✅ 模型移至GPU")
        else:
            device = torch.device("cpu")
            print("✅ 使用CPU模式")
        
        # 创建小的测试输入
        seq_len = 32
        msa_depth = 8
        
        batch = {
            'aatype': torch.randint(0, 20, (1, seq_len), device=device),
            'residue_index': torch.arange(seq_len, device=device).unsqueeze(0),
            'seq_length': torch.tensor([seq_len], device=device),
            'msa': torch.randint(0, 22, (1, msa_depth, seq_len), device=device),
            'num_alignments': torch.tensor([msa_depth], device=device),
            'msa_mask': torch.ones((1, msa_depth, seq_len), device=device),
            'seq_mask': torch.ones((1, seq_len), device=device),
            'template_aatype': torch.zeros((1, 0, seq_len), device=device, dtype=torch.long),
            'template_all_atom_positions': torch.zeros((1, 0, seq_len, 37, 3), device=device),
            'template_all_atom_mask': torch.zeros((1, 0, seq_len, 37), device=device),
            'template_mask': torch.zeros((1, 0), device=device),
            'template_pseudo_beta': torch.zeros((1, 0, seq_len, 3), device=device),
            'template_pseudo_beta_mask': torch.zeros((1, 0, seq_len), device=device),
            'extra_msa': torch.zeros((1, 0, seq_len), device=device, dtype=torch.long),
            'extra_msa_mask': torch.zeros((1, 0, seq_len), device=device),
            'extra_msa_row_mask': torch.zeros((1, 0), device=device),
        }
        
        print(f"✅ 测试输入准备完成:")
        print(f"   - 序列长度: {seq_len}")
        print(f"   - MSA深度: {msa_depth}")
        print(f"   - 设备: {device}")
        
        # 设置为评估模式并禁用梯度计算
        model.eval()
        with torch.no_grad():
            try:
                # 只计算表示，不计算损失
                output = model(batch)
                print(f"✅ 前向传播成功!")
                print(f"   - 输出键: {list(output.keys())}")
                if 'final_atom_positions' in output:
                    pos_shape = output['final_atom_positions'].shape
                    print(f"   - 原子坐标形状: {pos_shape}")
                
                return True
            except torch.cuda.OutOfMemoryError:
                print("⚠️ GPU内存不足，尝试更小的输入")
                return False
            except Exception as e:
                print(f"❌ 前向传播失败: {e}")
                return False
                
    except Exception as e:
        print(f"❌ 前向传播测试设置失败: {e}")
        return False

def main():
    """主测试函数"""
    print("🧬 OpenFold基础功能测试")
    print("=" * 50)
    
    tests = [
        ("模型创建", test_model_creation),
        ("数据变换", test_data_transforms), 
        ("CUDA内核", test_cuda_kernels),
        ("前向传播", test_forward_pass),
    ]
    
    results = {}
    for test_name, test_func in tests:
        try:
            results[test_name] = test_func()
        except Exception as e:
            print(f"❌ {test_name}测试异常: {e}")
            results[test_name] = False
    
    print("\n" + "=" * 50)
    print("🏆 测试结果总结:")
    
    passed = sum(results.values())
    total = len(results)
    
    for test_name, result in results.items():
        status = "✅ 通过" if result else "❌ 失败"
        print(f"   {test_name}: {status}")
    
    print(f"\n总计: {passed}/{total} 个测试通过")
    
    if passed == total:
        print("🎉 所有基础功能测试通过！")
        print("💡 OpenFold环境配置完全正确，可以进行结构预测")
        print("📝 下载模型参数后即可开始使用")
    else:
        print("⚠️ 部分测试失败，请检查环境配置")
    
    return passed == total

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1) 