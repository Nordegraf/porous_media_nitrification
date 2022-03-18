-- config for modelling contamination of a well
upwind_scheme = util.GetParam("--upwind", "full")

rhog = 998.23*9.81
numdays = 1000
tstop = numdays * 86400
tD = 86400 * 10
pump = -0.0002
entry = -0.00022

local well2D =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/well.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@Sandstone",
      type = "vanGenuchten",
      thetaS = 0.250, thetaR = 0.153,
      alpha = 0.79, n = 10.4,
      Ksat = 1.08},

    { uid = "@TouchetSiltLoam",
      type = "vanGenuchten",
      thetaS = 0.469, thetaR = 0.190,
      alpha = 0.50, n = 7.09,
      Ksat = 3.03},

    { uid = "@SiltLoam",
      type = "vanGenuchten",
      thetaS = 0.396, thetaR = 0.131,
      alpha = 0.423, n = 2.06,
      Ksat = 0.0496},
  },

  flow =
  {
    gravity = -9.81,            -- [ m s^{-2}], must be negative!
    density = 998.23,               -- [ kg m^{-3} ] saltwater density
    viscosity = 1.002e-3,           -- [ Pa s ]
    ammonium_diffusion = 1.86e-9,  -- [ m^2/s ]
    nitrate_diffusion = 1.7e-9,  -- [ m^2/s ]
    upwind = upwind_scheme, -- full, partial, no
  },
   medium =
   {
      {   subsets = {"Inner"},
          medium = "@SiltLoam",
      }
    },

  reactions =
  {
    rate = 0.6, --[1]
    molar_mass = 0.018039 --[kg/mol]
  },

  sources =
  {
    {cmp = "p", value = pump, subset = "Well", x = 1.0, y = 0.5},
    {cmp = "w_a", transport = pump, subset = "Well", x = 1.0, y = 0.5},
    {cmp = "w_n", transport = pump, subset = "Well", x = 1.0, y = 0.5},
  },

  initial =
  {
    { cmp = "w_a", value = 0.0 },
    { cmp = "w_n", value = 0.0 },
    { cmp = "p", value = "WellPressureStart" },
  },

  boundary =
  {
    -- Top
    {cmp = "p", type = "flux", bnd = "Top", inner="Inner", value = entry},
    {cmp = "w_a", type = "dirichlet", bnd = "Top", value = "manure"},

    -- Aquifer
    {cmp = "p", type = "dirichlet", bnd = "Aquifer", value = "WellPressureStart"},
    {cmp = "w_a", type = "dirichlet", bnd = "Aquifer", value = 0},
    {cmp = "w_n", type = "dirichlet", bnd = "Aquifer", value = 0},
  },

  linSolver =
  { type = "bicgstab",			-- linear solver type ["bicgstab", "cg", "linear"]
    precond =
    { type 		= "gmg",	                          -- preconditioner ["gmg", "ilu", "ilut", "jac", "gs", "sgs"]
      smoother 	= {type = "ilu", overlap = true},	-- gmg-smoother ["ilu", "ilut", "jac", "gs", "sgs"]
      cycle		= "V",		                          -- gmg-cycle ["V", "F", "W"]
      preSmooth	= 3,		                          -- number presmoothing steps
      postSmooth 	= 3,		                        -- number postsmoothing steps
      rap			= true,		                          -- comutes RAP-product instead of assembling if true
      baseLevel	= ARGS.numPreRefs,                -- gmg - baselevel
    },
    convCheck =
      { type		= "standard",
        iterations	= 30,		-- number of iterations
        absolute	= 0.5e-8,	-- absolut value of defact to be reached; usually 1e-8 - 1e-10 (must be stricter / less than in newton section)
        reduction	= 1e-7,		-- reduction factor of defect to be reached; usually 1e-7 - 1e-8 (must be stricter / less than in newton section)
        verbose		= true		-- print convergence rates if true
      }
  },

  time =
  {
    control	= "limex",
    start 	= 0.0,				      -- [s]  start time point
    stop	= tstop,			        -- [s]  end time point
    max_time_steps = 10000,		  -- [1]	maximum number of time steps
    dt		= 0.01,		          -- [s]  initial time step
    dtmin	= 0.001,	          -- [s]  minimal time step
    dtmax	= tstop/10,	            -- [s]  maximal time step
    dtred	= 0.5,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"w_a", "w_n", "p", "kr", "s", "q"},
  }

}

function manure(x, y, t)
  if t <= tD then
    return math.max(1.0-t/tD, 0.1)
  else
    return 0.1
  end
end

function rain(x, y, t)
  if t <= tD then
    return true, entry
  else
    return false, entry
  end
end

function WellPressureStart(x, y, t)
  return (1.0 - y) * rhog
end

return well2D