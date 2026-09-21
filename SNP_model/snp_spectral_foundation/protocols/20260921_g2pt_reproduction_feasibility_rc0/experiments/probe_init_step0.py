"""
B10 step-0 诊断探针：两种初始化下 HiGT block 的注意力统计与残差混合比。

只做 step 0，不训练。目的是判定 B10 的机制究竟落在哪一层——
是"注意力一开始就尖锐"（本探针已证伪），还是别的地方。

用法（先按 ENVIRONMENT_NOTES.zh-CN.md 建环境，并 clone 上游仓库）：
    cd <upstream G2PT clone>
    PYTHONPATH=. python probe_init_step0.py
"""
import math

import torch
import torch.nn as nn

from src.model.hierarchical_transformer.hierarchical_transformer import (
    HierarchicalTransformerUpdate,
)

D, H = 64, 4          # 论文配置：embedding dim 64, 4 heads
SEED = 0


def xavier_std(vocab, d=D):
    """xavier_uniform_ 对形状 (V, d): bound = sqrt(6/(d+V)), std = bound/sqrt(3)."""
    return math.sqrt(6.0 / (d + vocab)) / math.sqrt(3.0)


def make_block():
    torch.manual_seed(SEED)
    blk = HierarchicalTransformerUpdate(
        hidden=D, attn_heads=H, feed_forward_hidden=4 * D,
        norm=nn.LayerNorm(D), dropout=0.0, attention_dropout=0.0,
    )
    return blk.eval()


def _qk(blk, q_raw, k_raw):
    """复刻 forward 的路径：只有 q 过 norm_attn，k/v 原样进投影。"""
    dh = D // H
    B, Lq, _ = q_raw.shape
    Lk = k_raw.size(1)
    q_n = blk.norm_attn(q_raw)
    qh = blk.q_proj(q_n).view(B, Lq, H, dh).permute(0, 2, 1, 3)
    kh = blk.k_proj(k_raw).view(B, Lk, H, dh).permute(0, 2, 1, 3)
    vh = blk.v_proj(k_raw).view(B, Lk, H, dh).permute(0, 2, 1, 3)
    scores = (qh @ kh.transpose(-2, -1)) / math.sqrt(dh)
    return scores, vh


def attention_stats(blk, q_std, k_std, Lk, B=512, reps=8):
    ls, ent, mx = [], [], []
    for _ in range(reps):
        q_raw = torch.randn(B, 1, D) * q_std
        k_raw = torch.randn(B, Lk, D) * k_std
        with torch.no_grad():
            scores, _ = _qk(blk, q_raw, k_raw)
            p = torch.softmax(scores, dim=-1)
            ls.append(scores.std().item())
            ent.append((-(p * torch.log(p + 1e-12)).sum(-1)).mean().item())
            mx.append(p.max(-1).values.mean().item())
    mean = lambda v: sum(v) / len(v)
    return mean(ls), mean(ent), mean(mx)


def residual_ratio(blk, q_std, k_std, Lk, B=512, reps=8):
    """‖attn_output‖ / ‖q_raw‖，即 x = q + attn_output 里的混合比。"""
    out = []
    for _ in range(reps):
        q_raw = torch.randn(B, 1, D) * q_std
        k_raw = torch.randn(B, Lk, D) * k_std
        with torch.no_grad():
            scores, vh = _qk(blk, q_raw, k_raw)
            p = torch.softmax(scores, dim=-1)
            o = (p @ vh).permute(0, 2, 1, 3).reshape(B, 1, D)
            o = blk.o_proj(o)
            out.append((o.norm(dim=-1).mean() / q_raw.norm(dim=-1).mean()).item())
    return sum(out) / len(out)


def main():
    blk = make_block()
    # 论文 p=1e-5 档：5,510 SNP；代码里 snp_embedding 词表是 n_snps*3+2
    scales = {"snp": xavier_std(5510 * 3 + 2),
              "gene": xavier_std(253),
              "sys": xavier_std(20)}
    print("xavier 初始尺度（随词表大小分层）:")
    for k, v in scales.items():
        print(f"  {k:<5} std={v:.4f}")

    print("\n[1] step-0 注意力统计")
    print(f"{'Lk':>4} {'init':<9}{'logit std':>11}{'entropy':>10}{'ln(Lk)':>9}{'ratio':>8}{'max p':>8}")
    for Lk in (4, 8, 32):
        for name, (qs, ks) in (("default", (1.0, 1.0)),
                               ("xavier", (scales["gene"], scales["snp"]))):
            l, e, m = attention_stats(blk, qs, ks, Lk)
            print(f"{Lk:>4} {name:<9}{l:>11.4f}{e:>10.4f}"
                  f"{math.log(Lk):>9.4f}{e / math.log(Lk):>8.4f}{m:>8.4f}")

    print("\n[2] 残差混合比 ‖attn_out‖/‖q_raw‖")
    print(f"{'step':<28}{'default':>10}{'xavier':>10}{'ratio':>9}")
    for name, (qk_q, qk_k), Lk in (
        ("SNP->gene  (q=gene,k=SNP)", ("gene", "snp"), 8),
        ("gene->sys  (q=sys, k=gene)", ("sys", "gene"), 16),
        ("sys->gene  (q=gene,k=sys)", ("gene", "sys"), 6),
    ):
        rd = residual_ratio(blk, 1.0, 1.0, Lk)
        rx = residual_ratio(blk, scales[qk_q], scales[qk_k], Lk)
        print(f"{name:<28}{rd:>10.4f}{rx:>10.4f}{rd / rx:>9.1f}x")


if __name__ == "__main__":
    main()
