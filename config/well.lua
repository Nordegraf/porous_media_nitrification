-- config for modelling a drainage trench with constant groundwater flow

Trench2D_rho = 998.23
Trench2D_g = -9.81 -- must be negative!
rhog = (-1.0)*Trench2D_rho*Trench2D_g
numdays = 100
tstop = numdays * 86400 -- 100 days

local trench2D =
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
      alpha = 0.79/rhog, n = 10.4,
      Ksat = 1.08},

    { uid = "@TouchetSiltLoam",
      type = "vanGenuchten",
      thetaS = 0.469, thetaR = 0.190,
      alpha = 0.50/rhog, n = 7.09,
      Ksat = 3.03},

    { uid = "@SiltLoam",
      type = "vanGenuchten",
      thetaS = 0.396, thetaR = 0.131,
      alpha = 0.423/rhog, n = 2.06,
      Ksat = 0.0496},
  },

  flow =
  {
    gravity = Trench2D_g,            -- [ m s^{-2}], must be negative!
    density = 998.23,               -- [ kg m^{-3} ] saltwater density
    viscosity = 1.002e-3,           -- [ Pa s ]
    ammonia_diffusion = 18.8571e-6,  -- [ m^2/s ]
    nitrate_diffusion = 18.8571e-6,  -- [ m^2/s ]
  },
   medium =
   {
      {   subsets = {"Inner"},
          medium = "@SiltLoam",
      }
    },

  reactions =
  {
    nitrification = 0.2,
  },

  sources =
  {
    {cmp = "p", value = -0.0003, subset = "Well", x = 1.0, y = 0.5},
    --{cmp = "w_a", subset = "Well", x = 1.0, y = 0.5},
    --{cmp = "w_n", subset = "Well"}
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
    {cmp = "p", type = "flux", bnd = "Top", inner="Inner", value = -0.00009},
    {cmp = "w_n", type = "dirichlet", bnd = "Top", value = 1.0},

    -- Aquifer
    {cmp = "p", type = "dirichlet", bnd = "Aquifer", value = "WellAquiferBoundary" }
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
    dt		= 100,		          -- [s]  initial time step
    dtmin	= 0.001,	          -- [s]  minimal time step
    dtmax	= tstop/10,	            -- [s]  maximal time step
    dtred	= 0.5,			          -- [1]  reduction factor for time step
    tol 	= 1e-2,
  },

  output =
  {
    file = "./", -- must be a folder!
    data = {"w_a", "w_n", "p", "kr", "s", "q", "o"},
  }

}


function WellAquiferBoundary(x, y, t)
  return true, (1.0 - y) * rhog
end

function WellPressureStart(x, y, t)
  return (1.0 - y) * rhog
end

return trench2D
