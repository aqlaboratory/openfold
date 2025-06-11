#!/usr/bin/env python3
"""
OpenFold 模型参数下载脚本
自动下载所有必需的模型参数文件
"""

import urllib.request
import os
import sys
import time

def download_with_progress(url, target, description=""):
    """带进度条的下载函数"""
    def progress_hook(blocknum, blocksize, totalsize):
        readsofar = blocknum * blocksize
        if totalsize > 0:
            percent = readsofar * 100 / totalsize
            size_mb = totalsize / (1024*1024)
            read_mb = readsofar / (1024*1024)
            s = f"\r{description}: {percent:5.1f}% ({read_mb:6.1f}MB / {size_mb:6.1f}MB)"
            sys.stderr.write(s)
            if readsofar >= totalsize:
                sys.stderr.write("\n")
        else:
            read_mb = readsofar / (1024*1024)
            s = f"\r{description}: {read_mb:6.1f}MB 已下载"
            sys.stderr.write(s)

    try:
        print(f"开始下载: {url}")
        print(f"目标文件: {target}")
        os.makedirs(os.path.dirname(target), exist_ok=True)
        urllib.request.urlretrieve(url, target, progress_hook)
        print(f"✅ 下载完成: {target}")
        return True
    except Exception as e:
        print(f"❌ 下载失败: {e}")
        return False

def main():
    """主下载函数"""
    print("🧬 OpenFold 模型参数下载器")
    print("=" * 50)
    
    # 定义下载任务
    downloads = [
        {
            "url": "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar",
            "target": "openfold/resources/params/alphafold_params_2022-12-06.tar",
            "description": "AlphaFold2参数",
            "size": "4.7GB"
        },
        {
            "url": "https://huggingface.co/nz/OpenFold/resolve/main/finetuning_1.pt",
            "target": "openfold/resources/openfold_params/finetuning_1.pt", 
            "description": "OpenFold参数-1",
            "size": "~200MB"
        },
        {
            "url": "https://huggingface.co/nz/OpenFold/resolve/main/finetuning_2.pt",
            "target": "openfold/resources/openfold_params/finetuning_2.pt",
            "description": "OpenFold参数-2", 
            "size": "~200MB"
        },
        {
            "url": "https://huggingface.co/nz/OpenFold/resolve/main/finetuning_3.pt",
            "target": "openfold/resources/openfold_params/finetuning_3.pt",
            "description": "OpenFold参数-3",
            "size": "~200MB"
        },
        {
            "url": "https://huggingface.co/nz/OpenFold/resolve/main/finetuning_4.pt", 
            "target": "openfold/resources/openfold_params/finetuning_4.pt",
            "description": "OpenFold参数-4",
            "size": "~200MB"
        },
        {
            "url": "https://huggingface.co/nz/OpenFold/resolve/main/finetuning_5.pt",
            "target": "openfold/resources/openfold_params/finetuning_5.pt",
            "description": "OpenFold参数-5",
            "size": "~200MB"
        }
    ]
    
    # 显示下载计划
    print("📋 下载计划:")
    total_files = len(downloads)
    for i, item in enumerate(downloads, 1):
        print(f"   {i}. {item['description']} ({item['size']})")
    print()
    
    # 检查已存在的文件
    print("🔍 检查已存在的文件...")
    to_download = []
    for item in downloads:
        if os.path.exists(item['target']):
            size_mb = os.path.getsize(item['target']) / (1024*1024)
            print(f"   ✅ {item['description']}: 已存在 ({size_mb:.1f}MB)")
        else:
            to_download.append(item)
            print(f"   ⏳ {item['description']}: 需要下载")
    
    if not to_download:
        print("\n🎉 所有文件已存在，无需下载！")
        return True
    
    print(f"\n📥 需要下载 {len(to_download)} 个文件")
    
    # 开始下载
    success_count = 0
    for i, item in enumerate(to_download, 1):
        print(f"\n[{i}/{len(to_download)}] {item['description']}")
        success = download_with_progress(
            item['url'], 
            item['target'], 
            item['description']
        )
        if success:
            success_count += 1
        time.sleep(1)  # 避免请求过于频繁
    
    # 总结
    print("\n" + "=" * 50)
    print("📊 下载总结:")
    print(f"   成功: {success_count}/{len(to_download)}")
    print(f"   失败: {len(to_download) - success_count}")
    
    if success_count == len(to_download):
        print("\n🎉 所有参数文件下载完成！")
        print("📝 下一步: 解压AlphaFold参数")
        print("   cd openfold/resources/params")
        print("   tar -xf alphafold_params_2022-12-06.tar")
        return True
    else:
        print("\n⚠️ 部分文件下载失败")
        print("💡 建议:")
        print("   1. 检查网络连接")
        print("   2. 使用代理或镜像")
        print("   3. 手动下载失败的文件")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 