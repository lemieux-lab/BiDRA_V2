(pwd() != @__DIR__) && cd(@__DIR__) # allow starting app from bin/ dir

using BiDRADashboardV2
const UserApp = BiDRADashboardV2
BiDRADashboardV2.main()
