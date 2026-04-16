import argparse
import os
import math

def compare_text(file0,file1):
    with open(file0) as f0:
        with open(file1) as f1:
            for l0, l1 in zip(f0.readlines(),f1.readlines()):
                if l0!=l1:
                    print(l0)
                    print(l1)
                    assert(False)
                    return
    assert(True)

def compare_nearly(file0,file1,thresh):
    with open(file0) as f0:
        with open(file1) as f1:
            for l0, l1 in zip(f0.readlines(),f1.readlines()):
                # テキストで一致なら pass
                if l0!=l1:
                    # テキスト一致しない場合は数値比較
                    v0 = map(float,l0.strip().split())
                    v1 = map(float,l1.strip().split())
                    for u0,u1 in zip(v0,v1):
                        if math.fabs(u0-u1) > thresh:
                            print(l0)
                            print(l1)
                            assert(False)
                            return
    assert(True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='FastBN Check Program')
    parser.add_argument('result_file', help='result file')
    parser.add_argument('reference_file', help='reference file')
    parser.add_argument('--type', '-t', default='text', help='compare by text or nearly')
    args = parser.parse_args()
    if args.type=="text":
        print(f"compare {args.result_file} with {args.reference_file} by text")
        compare_text(args.result_file,args.reference_file)
    if args.type=="nearly":
        print(f"compare {args.result_file} with {args.reference_file} by number")
        compare_nearly(args.result_file,args.reference_file,1.0e-6)
