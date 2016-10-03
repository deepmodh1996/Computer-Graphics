(Material matRED
	(color 0.75 0.25 0.25 )
	(ma 0.5 ) (ml 0.6 ) (ms 0.1 ) (msp 150 ) (msm 1 ) (mt 0 ) (mtn 1 ))
(Material matBLUE
	(color 0.25 0.25 0.75 )
	(ma 0.5 ) (ml 0.6 ) (ms 0.1 ) (msp 20 ) (msm 1 ) (mt 0 ) (mtn 1 ))
(Material matGREEN
	(color 0.25 0.75 0.25 )
	(ma 0.5 ) (ml 0.6 ) (ms 0.1 ) (msp 150 ) (msm 1 ) (mt 0 ) (mtn 1 ))
(Material matWHITE
	(color 0.75 0.75 0.75 )
	(ma 0.5 ) (ml 0.6 ) (ms 0.1 ) (msp 40 ) (msm 1 ) (mt 0 ) (mtn 1 ))
(Material matGlass
	(color 0.9 0.9 0.9 )
	(ma 0.0 ) (ml 0.2 ) (ms 0.1 ) (msp 40 ) (msm 1 ) (mt 0.9 ) (mtn 1.007))
(Light amb (type 0 ) (color 0.2 0.2 0.2 ) )
(Light pt1 (type 2 ) (color 20.8 20.8 20.8 )  (falloff 1 )  (deaddistance 0.1 ) )
(Camera cam
	(perspective 1 )
	(l -1 ) (r 1 ) (b -1 ) (t 1 ) (n 3 )) #no far plane
(Sphere spTOP (radius 100000.0 ) (material matGREEN ))
(Sphere spWALL (radius 100000.0 ) (material matWHITE ))
(Sphere spLEFT (radius 100000.0 ) (material matBLUE ))
(Sphere spRIGHT (radius 100000.0 ) (material matRED ))
(Sphere spGlass (radius 1 ) (material matGlass ))
(Cylinder clGlass (radius 1 ) (height 4) (material matGlass ))
(Sphere spNEW (radius 1 ) (material matGREEN ))
(G awesome
	(I cam1 cam )
	(I dir1 dir (Xform 0 0 -0.5 0
		0 0 -0.5 0
		0 0 0.5 0
		0 0 0 0 ) )
	(I amb1 amb )
	(I pt1_I pt1 (T 2 2 -45 ) )
	(I spGlass_I spNEW (T 4 0 -50 ) )
	(I cl2 clGlass (R 90 1 0 0) (R 30 0 0 1) (T -0.60 1.4 -10 ) )
	(I spB_I spWALL (T 0 -100009.0 -3 ) )
	(I spT_I spTOP (T 0 100009.0 -3 ) )
	(I spF_I spWALL (T 0 0 -100100.0 ) )
	(I spN_I spWALL (T 0 0 100000.0 ) )
	(I spL_I spLEFT (T -100009.0 0 -3 ) )
	(I spR_I spRIGHT (T 100009.0 0 -3 ) )
)
(Render awesome )
