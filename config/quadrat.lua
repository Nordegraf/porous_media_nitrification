rho = 998.23
g = -9.81 -- must be negative!
rhog = (-1.0)*rho*g
tstop = 10

local quadrat =
{
  -- The domain specific setup
  domain =
  {
    dim = 2,
    grid = "grids/quadrat.ugx",
    numRefs = ARGS.numRefs,
    numPreRefs = ARGS.numPreRefs,
  },

  -- medium parameters for vanGenuchten Model
  parameter = {
    { uid = "@Sand",
      type = "vanGenuchten",
      thetaS = 0.37, thetaR = 0.043,
      alpha = 0.087/rhog, n = 1.58,
      Ksat = 1}
  },

  flow =
  {
    gravity = g,            -- [ m s^{-2}], must be negative!
    density = 998.23,               -- [ kg m^{-3} ] saltwater density
    viscosity = 1.002e-3,           -- [ Pa s ]
    ammonia_diffusion = 18.8571e-6,  -- [ m^2/s ]
    nitrate_diffusion = 18.8571e-6,  -- [ m^2/s ]
    nitrite_diffusion = 18.8571e-6  -- [ m^2/s ]
  },

  reactions =
  {
    ammonia_rate = 0.2,
    nitrite_rate = 0.1
  },

  medium =
  {
    { subsets = {"Inner"},
      medium = "@Sand"
    },
  },

  initial =
  {
    { cmp = "c_am", value = 0.0 },
    { cmp = "c_it", value = 0.0 },
    { cmp = "c_at", value = 0.0 },
    { cmp = "p", value = "PressureStart" },
  },

  boundary =
  {
    { cmp = "p", type = "flux", bnd = "Upper", inner="Inner", value = -1},
    { cmp = "p", type = "dirichlet", bnd = "Outflow", value = "PressureStart"},
    { cmp = "c_am", type = "dirichlet", bnd = "Upper", value = 1.0}
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
    max_time_steps = 1000,		  -- [1]	maximum number of time steps
    dt		= 1000,		          -- [s]  initial time step
    dtmin	= 0.001,	          -- [s]  minimal time step
    dtmax	= tstop/10,	            -- [s]  maximal time step
    dtred	= 0.5,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"c_am", "c_it", "c_at", "p", "kr", "s", "q", "pc", "o"},
  },
}

function PressureStart(x, y, t)
  return (1.0 - y) * rhog
end

return quadrat