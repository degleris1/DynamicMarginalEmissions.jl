abstract type AbstractDevice end

not_implemented() = error("Must implement this!")


# Device Data
get_cost(d::AbstractDevice, p) = not_implemented()

make_aux_vars(d::AbstractDevice) = not_implemented()

make_constraints(d::AbstractDevice, p, aux_vars) = make_constraints(d, p)
make_constraints(d::AbstractDevice, p) = not_implemented()


# Dimensions, indices
get_time_horizon(d::AbstractDevice) = not_implemented()

get_num_aux(d::AbstractDevice) = 0

get_dims(d::AbstractDevice) = not_implemented()


# KKTs
get_device_dL_dx(d::AbstractDevice, p, duals, aux_vars) = get_device_dL_dx(d, p, duals)
get_device_dL_dx(d::AbstractDevice, p, duals) = not_implemented()

get_device_dL_daux(d::AbstractDevice, p, duals, aux_vars) = not_implemented()
get_device_dL_daux(d::AbstractDevice, p, duals, aux_vars::Nothing) = Real[]

get_device_comp_slack(d::AbstractDevice, p, duals, aux_vars) = get_device_comp_slack(d, p, duals)
get_device_comp_slack(d::AbstractDevice, p, duals) = not_implemented()


# Jacobian
jacobian_kkt_z_device(d::AbstractDevice, p, duals, aux_vars) = jacobian_kkt_z_device(d::AbstractDevice, p, duals)
jacobian_kkt_z_device(d::AbstractDevice, p, duals) = not_implemented()