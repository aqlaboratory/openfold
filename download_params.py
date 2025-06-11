#!/usr/bin/env python3
import urllib.request
import os
import sys

def download_with_progress(url, target):
    """下载文件并显示进度"""
    def progress_hook(blocknum, blocksize, totalsize):
        readsofar = blocknum * blocksize
        if totalsize > 0:
            percent = readsofar * 1e2 / totalsize
            s = f"\r{percent:5.1f}% {readsofar / (1024*1024):8.1f} MB / {totalsize / (1024*1024):8.1f} MB"
            sys.stderr.write(s)
            if readsofar >= totalsize:
                sys.stderr.write("\n")
        else:
            sys.stderr.write(f"\r已下载 {readsofar / (1024*1024):8.1f} MB")

    print(f'开始下载: {url}')
    print(f'目标文件: {target}')
    
    os.makedirs(os.path.dirname(target), exist_ok=True)
    
    try:
        urllib.request.urlretrieve(url, target, progress_hook)
        print(f'\n下载完成: {target}')
        print(f'文件大小: {os.path.getsize(target) / (1024*1024*1024):.2f} GB')
        return True
    except Exception as e:
        print(f'\n下载失败: {e}')
        return False

if __name__ == "__main__":
    url = 'https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar'
    target = 'openfold/resources/params/alphafold_params_2022-12-06.tar'
    
    success = download_with_progress(url, target)
    if success:
        print("AlphaFold参数下载成功！")
    else:
        print("AlphaFold参数下载失败！")
        sys.exit(1) 