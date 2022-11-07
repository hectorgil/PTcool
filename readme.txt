#Hector Gil Marin 7 Nov 22. PTcool readme file
#Terms of the perturbation theory file produced by PTcool

#0 k (in h/Mpc)

#1 Plin (in [Mpc/h]^3). Definition as <delta_k' delta_k > = (2pi)^3 delta^D(k+k') * P(k)
#2 P22_delta_delta (Eq. A5 of 1209.3771)
#3 P22_delta_theta (like Eq. A5 of 1209.3771, but changing of the F2 squared by G2)
#4 P22_theta_theta (like Eq. A5 of 1209.3771, but changing F2 by G2)
#5 P33_delta_delta (second term of Eq. A8 of 1209.3771)
#6 P33_delta_theta (like second term of Eq. A8 of 1209.3771, but chaning one of the F3 squared by G3)
#7 P33_theta_theta (like second term of Eq. A8 of 1209.3771, but chaning F3 by G3)
#8 2*P13_delta (Eq. A4 of 1209.3771)
#9 2*P13_theta (like Eq. A4 of 1209.3771, but changing F3 by G3)
#10 2*P15_delta (Eq. A6 of 1209.3771)
#11 2*P15_theta (like Eq. A6 of 1209.3771, but changing F5 by G5)
#12 2*P24_delta_delta (Eq A7 of 1209.3771)
#13 2*P24_delta_delta (like Eq A7 of 1209.3771, but changing F2 and F4 by G2 and G4)
#14 2*P24_delta_theta_1 (hybrid form 1 of P24 cross term: F2*F4 -> G2*F4 in Eq. A7 of 1209.3771)
#15 2*P24_delta_theta_2 (hybrid form 2 of P24 cross term: F2*F4 -> F2*G4 in Eq. A7 of 1209.3771) 

	#Real Space P_delta_delta, P_delta_theta and P_theta_theta terms
	Note that for SPT 1L P_delta_delta = Plin + (P22_delta_delta+2*P13_delta) (Eq. A2 of 1209.3771)
	Note that for SPT 2L P_delta_delta = Plin + (P22_delta_delta+2*P13_delta) + (2*P15_delta_delta+2*P24_delta_delta+P33_delta_delta + P13_delta_delta*P13_delta_delta/(4*Plin) ) (Eq. A3 of 1209.3771)
	Similarly P_theta_theta can be found by changing all P??_delta_delta terms by P??_theta_theta.
	For the cross delta-theta,
	SPT 1L P_delta_theta = Plin + (P22_delta_theta+0.5*(2*P13_delta+2*P13_theta))
	SPT 2L P_delta_theta = Plin + (P22_delta_theta+0.5*(2*P13_delta+2*P13_theta)) + (P33_delta_theta + 0.5*(2*P24_delta_theta_1+2*P24_delta_theta_2) + 0.5*(2*P15_delta+2*P15_theta) + 0.25*(2*P13_delta+2*P13_theta)*(2*P13_delta+2*P13_theta)/(4*Plin))

	For the 1L RPT, P_delta_delta = ( Plin + P22_delta_theta )*exp(2*P13_delta/Plin) (Eq. B27 of 1209.3771)
	For the 2L RPT, P_delta_delta = ( Plin + P22_delta_theta + P33_delta_delta)*N2^2 (Eq. B38 of 1209.3771)
	where N2 = cosh( sqrt(2*P15_delta/Plin) ) + 0.5*2*P13_delta/Plin*sqrt(Plin/(2*P15_delta))*sinh(sqrt(2*P15_delta/Plin))  (Eq. B39 of 1209.3771)
	Similarly P_theta_theta can be found by changing all P??_delta_delta terms by P??_theta_theta.
	For the cross delta-theta,
	1L RPT, P_delta_delta P_delta_theta = (Plin + P22_delta_theta)*exp(0.5*(2*P13_delta+2*P13_theta)/Plin)
	2L RPT, P_delta_delta P_delta_theta = (Plin + P22_delta_theta + P33_delta_theta)*N2_cross
	where N2_cross is like N2 but changing 2*P13_delta and 2*P15_theta terms by 0.5*(2*P13_delta+2*P13_theta) and 0.5*(2*P15_delta+2*P15_theta), respectively.

	(see Eq. B1 and B9 of 1407.5668 on how to add the following terms)

#16 Pb2_delta (Eq. B2 of 1407.5668)
#17 Pb2_theta (Eq. B10 of 1407.5668)
#18 Pbs2_delta (Eq. B3 of 1407.5668)
#19 Pbs2_theta (Eq. B11 of 1407.5668)
#20 Pb2s2 (Eq. B4 of 1407.5668)
#21 Pbs22 (Eq. B5 of 1407.5668) 
#22 Pb22 (Eq. B6 of 1407.5668)
#23 Plin*sigma3^2 (see Eq. B7 of 1407.5668)

	(see Appendix A of 1006.0699 for the definition of A and B following terms)

#24 A11
#25 A12
#26 A22
#27 A23
#28 A33
#29 B1_11
#30 B1_12
#31 B1_21
#32 B1_22
#33 B2_11
#34 B2_12
#35 B2_21
#36 B2_22
#37 B3_12
#38 B3_21
#39 B3_22
#40 B4_22

	(see Eq B15 and B16 of 1407.5668)
	A_TNS(k,mu,b1) = A11*b1*b1*f*mu^2 + A12*b1*f*f*mu^2 + A22*b1*f*f*mup^4 + A23*f*f*f*mu^4 + A33*f*f*f*mu^6
	B_TNS(k,mu,b1) = B1_11*b1*b1*f*f*mu^2 + B1_12*b1*f*f*f*mu^2 + B1_21*b1*f*f*f*mu^2 + B1_22*f*f*f*f*mu^2 + B2_11*b1*b1*f*f*mu^4 + B2_12*b1*f*f*f*mu^4 + B2_21*b1*f*f*f*mu^4 + B2_22*f*f*f*f*mu^4 + B3_12*b1*f*f*f*mu^6 + B3_21*b1*f*f*f*mu^6 + B3_22*f*f*f*f*mu^6 + B4_22*f*f*f*f*mu^8

#41 sigma8 (value of sigma8 computed from Plin)

Note you need to appropiatelly to re-scale each of these terms above if the amplitude of the linear power spectrum changes. Plin~sigma8^2, P13 and P22 ~sigma8^4, etc
