module Constants
@enum SPIN unpolarized=1 polarized=2
@enum RELATIVISTIC relativistic=1 non_relativistic=0
@enum KIND exchange=0 correlation=1 excorr=2 kinetic=3
""" Rung on Jacob's ladder """
@enum FAMILY unknown=-1 lda=1 gga=2 mgga=4 lca=8 oep=16 hybrid_gga=32 hybrid_mgga=64
@enum(FLAGS, exc=(1<<0), vxc=(1<<1), fxc=(1<<2), kxc=(1<<3), lxc=(1<<4),
      D1=(1<<5), D2=(1<<6), D3=(1<<7),
      hyb_cam=(1<<8), hyb_camy=(1<<9), VV10=(1<<10), hyb_lc=(1<<11), hyb_lcy=(1<<12),
      stable=(1<<13), development=(1<<14))
@enum TAU explicit=0 expansion=1
end
