abstract type AbstractDevice end

not_implemented() = error("Must implement this!")

# Device Data
get_cost(d::AbstractDevice, p) = not_implemented()
get_constraints(d::AbstractDevice, p) = not_implemented()

# Dimensions, indices
get_time_horizon(d::AbstractDevice) = not_implemented()
get_dims(d::AbstractDevice) = not_implemented()

# KKTs
get_device_dL_dx(d::AbstractDevice, p, duals) = not_implemented()
get_device_comp_slack(d::AbstractDevice, p, duals) = not_implemented()

# Jacobian
jacobian_kkt_z_device(d::AbstractDevice, p, duals) = not_implemented()