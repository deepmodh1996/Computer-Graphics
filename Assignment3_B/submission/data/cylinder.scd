(Material matRED
	(color 0.75 0.25 0.25 )
	(ma 0.5 ) (ml 0.6 ) (ms 0.1 ) (msp 150 ) (msm 1 ) (mt 0 ) (mtn 1 ))
(Light amb (type 0 ) (color 0.2 0.2 0.2 ) )
(Light pt1 (type 2 ) (color 20.8 20.8 20.8 )  (falloff 1 )  (deaddistance 0.1 ) )
(Camera cam
	(perspective 1 )
	(l -1 ) (r 1 ) (b -1 ) (t 1 ) (n 3 )) #no far plane
(Cylinder cyl (radius 1 ) (height 1.5) (material matRED ))
(G awesome
	(I cam1 cam )
	(I amb1 amb )
	(I pt1_I pt1 (T 2 2 -5 ) )
	(I cyl1 cyl (S 1 0.5 1) (R 60 1 1 0) (T 0 1 -10 ) )
)
(Render awesome )
