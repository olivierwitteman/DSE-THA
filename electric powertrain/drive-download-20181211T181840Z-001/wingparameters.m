wingarea = input("What's the wing area? ")
span = sqrt(wingarea*a.AR)
spantakenupbypropulsors = p.b_dp*span
propulsordiameter = span*p.b_dp/(p.N2*(1+p.dy))