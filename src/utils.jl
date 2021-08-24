function make_dynamic(net::PowerNetwork, T, P, C, dyn_gmax, η_c, η_d)
	fqs = [net.fq for t in 1:T]
	fls = [net.fl for t in 1:T]
	pmaxs = [net.pmax for t in 1:T]
	return DynamicPowerNetwork(fqs, fls, pmaxs, dyn_gmax, net.A, net.B, P, C, T; η_c=η_c, η_d=η_d)
end