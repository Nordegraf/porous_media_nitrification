-- author: Niklas Conen
ug_load_script("ug_util.lua")

-- namespace
util.unsat = util.unsat or {}

local json = require("json")

-- object oriented clustering of the Approximation Space, Element Discretisation and
-- Domain Discretisation
ProblemDisc = {}

function ProblemDisc:new(problemDesc, dom)
    assert(problemDesc, "No Input defined")
    problemObj = {}
    setmetatable(problemObj, self)
    self.__index = self

    self.problem = problemDesc
    self.domain = dom
    self.cmp = problemDesc.flow.cmp

    self.gravity = ConstUserVector(0.0)
    self.gravity:set_entry(problemDesc.domain.dim-1, problemDesc.flow.gravity)

    self.modelMap = nil
    if problemDesc.parameter ~= nil then
        self.modelMap = ProblemDisc:CreateModelMap(problemDesc.parameter)
    end

    return problemObj
end

function ProblemDisc:CreateApproxSpace()
    -- documentation for ug function ApproximationSpace: line 277 in /ug4/ugcore/ugbase/lib_disc/function_spaces/approximation_space.h
    self.approxSpace = ApproximationSpace(self.domain)
    self.approxSpace:add_fct("p", "Lagrange", 1)
    self.approxSpace:add_fct("w_a", "Lagrange", 1)
    self.approxSpace:add_fct("w_n", "Lagrange", 1)
    self.approxSpace:init_levels()
    self.approxSpace:init_top_surface()
    self.approxSpace:print_statistic()
    return self.approxSpace
end

-- Element discretisation
function ProblemDisc:CreateElemDisc(subdom, medium)

    -- Creates the elememt discretisation for a given medium
	local elemDisc = {}
    -- flow equation
	elemDisc["p"] = ConvectionDiffusion("p", subdom, "fv1")
	-- transport equations
    elemDisc["w_a"] = ConvectionDiffusion("w_a", subdom, "fv1")
    elemDisc["w_n"] = ConvectionDiffusion("w_n", subdom, "fv1")

    -- constant viscosity and density
    viscosity = self.problem.flow.viscosity
    density = self.problem.flow.density

    local porosity = nil    -- phi
    for i, param in ipairs(self.problem.parameter) do
        if param.uid == medium.porosity then
            porosity = param.thetaS     -- the porosity is equal to the saturated water content
        end
    end

    Ksat = nil
    porosity = nil
    for i, param in ipairs(self.problem.parameter) do
        if param.uid == medium.medium then
            Ksat = param.Ksat
            porosity = param.thetaS
        end
    end

    local capillary = -1.0*elemDisc["p"]:value()

    -- used for switching off the van Genuchten Model for debugging
    local conductivity = nil
    local saturation = nil
    if true then
        -- the vanGenuchten model is calculated using the Richards Plugin
        conductivity = ProblemDisc:conductivity(medium.medium) -- k(S)
        saturation = ProblemDisc:saturation(medium.medium) -- S
        conductivity:set_capillary(capillary)
        saturation:set_capillary(capillary)
    else
        conductivity = Ksat
        saturation = 1.0
    end

    local r_scale = -86400 -- scaling for time unit in K_sat
    local permeability = ScaleAddLinkerMatrix()
    -- needs scaling depending on hydraulic conductivities scale in the Richardsplugin
    permeability:add(conductivity, self.problem.flow.viscosity/(r_scale*self.problem.flow.density*self.problem.flow.gravity))
    print("saturated permeability: "..(Ksat*self.problem.flow.viscosity/(r_scale*self.problem.flow.density*self.problem.flow.gravity)))

    -- Darcy Velocity
    -- $\vec q := -k*k(p)/mu (\grad p - \rho \vec g)$
    local DarcyVelocity = DarcyVelocityLinker()
    DarcyVelocity:set_viscosity(viscosity)
    DarcyVelocity:set_permeability(permeability)
    DarcyVelocity:set_pressure_gradient(elemDisc["p"]:gradient())
    DarcyVelocity:set_density(density)
    DarcyVelocity:set_gravity(self.gravity)
    self.DarcyVelocity = DarcyVelocity

    local volufrac = ScaleAddLinkerNumber()
    volufrac:add(porosity, saturation)

	-----------------------------------------
	-- Equation [1] - Fluid Flow Equation
	-----------------------------------------
	-- $\partial_t (\Phi S_w \rho_w)
	--		+ \nabla \cdot (\rho_w \vec{v}_w) = \rho_w \Gamma_w$

	-- fluid storage: \Phi S_w \rho_w
    local storage = ScaleAddLinkerNumber()
    storage:add(volufrac, density)
	-- flux of the fluid phase \rho_w \vec{v}_w
    local fluidFlux = ScaleAddLinkerVector()
    fluidFlux:add(density, DarcyVelocity)

    elemDisc["p"]:set_mass(storage)
    elemDisc["p"]:set_flux(fluidFlux)
    elemDisc["p"]:set_mass_scale(0.0)

    -----------------------------------------
	-- Equation [2] - Nitrate transport
	-----------------------------------------
    -- \frac{\partial}{\partial t}(\phi S_w \rho \omega_{A}) + \nabla  (\bvec{q} \rho \omega_{A} - \phi S_w \rho D \nabla \omega_{A}) = \phi S_w \rho \Gamma - \phi S_w \rho k

    local nitrate_diffusion = ScaleAddLinkerMatrix()
    nitrate_diffusion:add(volufrac, self.problem.flow.nitrate_diffusion*self.problem.flow.density)

    elemDisc["w_n"]:set_mass_scale(storage)
    elemDisc["w_n"]:set_velocity(fluidFlux)
    elemDisc["w_n"]:set_diffusion(nitrate_diffusion)

    -- -phi * S * M/rho * k * omega^A
    function nitrification_switch(p)
        if p <= 0 then
            return 1
        else
            return 0
        end
    end

    function dnitrification_switch(p)
        return 0
    end

    local heaviside_nit = LuaUserFunctionNumber("nitrification_switch", 1)
    heaviside_nit:set_deriv(0, "dnitrification_switch")
    heaviside_nit:set_input(0, elemDisc["p"]:value())

    -- phi * S * M/rho * k * omega^A
    local nitrification = ScaleAddLinkerNumber()
    nitrification:add(self.problem.reactions.rate * self.problem.reactions.molar_mass / self.problem.flow.density, volufrac*elemDisc["w_a"]:value()*heaviside_nit)

    elemDisc["w_n"]:set_source(nitrification)

	-----------------------------------------
	-- Equation [3] - Ammonium transport
	-----------------------------------------
    -- \frac{\partial}{\partial t}(\phi S_w \rho \omega_{N}) + \nabla  (\bvec{q} \rho \omega_{N} - \phi S_w \rho D \nabla \omega_{N}) = \phi S_w \rho \Gamma + \phi S_w \rho k

    local am_diffusion = ScaleAddLinkerMatrix()
    am_diffusion:add(volufrac, self.problem.flow.ammonium_diffusion*self.problem.flow.density)

    elemDisc["w_a"]:set_mass_scale(storage)
    elemDisc["w_a"]:set_velocity(fluidFlux)
    elemDisc["w_a"]:set_diffusion(am_diffusion)

    -- - phi * S * M/rho * k * omega^A
    local ammonium_oxidation = ScaleAddLinkerNumber()
    ammonium_oxidation:add(-1.0, nitrification)
    elemDisc["w_a"]:set_source(ammonium_oxidation)

    if self.problem.flow.upwind == "full" then
        elemDisc["w_a"]:set_upwind(FullUpwind())
        elemDisc["w_n"]:set_upwind(FullUpwind())
    elseif self.problem.flow.upwind == "partial" then
        elemDisc["w_a"]:set_upwind(PartialUpwind())
        elemDisc["w_n"]:set_upwind(PartialUpwind())
    end

    local si = self.domain:subset_handler():get_subset_index(subdom)

    self.CompositeCapillary:add(si, capillary)
    self.CompositeConductivity:add(si, conductivity)
    self.CompositeSaturation:add(si, saturation)
    self.CompositeDarcyVelocity:add(si, DarcyVelocity)

    print("Created Element Discretisation for Subset ", subdom)

    return elemDisc
end


function ProblemDisc:CreateDomainDisc(approxSpace)
    self.domainDisc = DomainDiscretization(approxSpace)

    self.CompositeCapillary = CompositeUserNumber(true)
    self.CompositeConductivity = CompositeUserNumber(false)
    self.CompositeSaturation = CompositeUserNumber(false)
    self.CompositeDarcyVelocity = CompositeUserVector(false)

    for i,medium in ipairs(self.problem.medium) do
        local elemDisc = nil

        for j, subset in ipairs(medium.subsets) do
            elemDisc = self:CreateElemDisc(subset, medium)

            self.domainDisc:add(elemDisc["p"])
            self.domainDisc:add(elemDisc["w_a"])
            self.domainDisc:add(elemDisc["w_n"])
        end
    end

    -- Create Boundary Conditions
    local dirichletBnd = nil
    local neumannBnd = {}
    for i, v in ipairs(self.problem.boundary) do

        if v.type == "dirichlet" then
            -- Dirichtlet boundary
            dirichletBnd = dirichletBnd or DirichletBoundary()
            dirichletBnd:add(v.value, v.cmp, v.bnd)
            print("Added Dirichlet Boundary with value " .. v.value .. " for " .. v.cmp .. " on subset " .. v.bnd)
        end

        if v.type == "flux" then
            -- Neumann-type
            neumannBnd[v.cmp] = neumannBnd[v.cmp] or NeumannBoundary(v.cmp, "fv1")
            neumannBnd[v.cmp]:add(v.value, v.bnd, v.inner)
            print("Added Neumann Boundary with value " .. v.value .. " for " .. v.cmp .. " on subset " .. v.bnd)
        end

    end

    if (dirichletBnd) then  self.domainDisc:add(dirichletBnd) end
    if (neumannBnd["p"]) then  self.domainDisc:add(neumannBnd["p"]) end
    if (neumannBnd["w_a"]) then  self.domainDisc:add(neumannBnd["w_a"]) end
    if (neumannBnd["w_n"]) then  self.domainDisc:add(neumannBnd["w_n"]) end

    -- Adding Dirac Sources
    if self.problem.sources ~= nil then
        for i, v in ipairs(self.problem.sources) do

            local source = DiracSourceDisc(v.cmp, v.subset)

            local location = nil
            if self.problem.domain.dim == 2 and v.x and v.y then
                location = Vec2d(v.x, v.y)
            elseif self.problem.domain.dim == 3 and v.x and v.y and v.z then
                location = Vec3d(v.x, v.y, v.z)
            else
                print("source coordinates not given")
            end

            if v.value then source:add_source(v.value, location)
            elseif v.transport then source:add_transport_sink(v.transport)
            end
            self.domainDisc:add(source)

            if v.value ~= nil then
                print("Added DiracPointSource with value " .. v.value .. " on subset " .. v.subset)
            elseif v.transport ~= nil then
                print("Added DiracPointSource with transport value" ..v.transport.. "on subset " .. v.subset)
            end
        end
    end

    print("Created Domain Discretisation")
    return self.domainDisc
end


function ProblemDisc:CreateVTKOutput()
    -- generating the vtk output
    -- for future reference; util.Balance()
    for i, v in ipairs(self.problem.output.data) do
        -- concentration and pressure
        if v == "p" then
           self.vtk:select_nodal(GridFunctionNumberData(self.u, v), v)
        -- concentrations
        elseif v == "w_a" then
            self.vtk:select_nodal(GridFunctionNumberData(self.u, v), v)
        elseif v == "w_n" then
            self.vtk:select_nodal(GridFunctionNumberData(self.u, v), v)
        -- relative conductivity
        elseif v == "kr" then
            self.vtk:select_element(self.CompositeConductivity, v)
        -- saturation
        elseif v == "s" then
            self.vtk:select_element(self.CompositeSaturation, v)
        -- darcy velocity
        elseif v == "q" then
            self.vtk:select(self.CompositeDarcyVelocity, v)
        -- capillary pressure
        elseif v == "pc" then
            self.vtk:select(self.CompositeCapillary, v)
        end
    end

end

function ProblemDisc:CreateModelMap(paramDesc)
    local modelMap = {}
    for i, medium in ipairs(paramDesc) do
        if medium.type == "vanGenuchten" then
            medium.alpha= medium.alpha / (-1.0 * self.problem.flow.gravity * self.problem.flow.density)
            modelMap[medium.uid] = CreateVanGenuchtenModel(json.encode(medium))
        elseif medium.type == "const" then
            modelMap[medium.uid] = medium.value
        end
    end
    return modelMap
end


function ProblemDisc:SetInitialData()
    for i, initial in ipairs(self.problem.initial) do
        print("Added Initial Value " .. initial.cmp .. " = " .. initial.value)
        Interpolate(initial.value, self.u, initial.cmp)
    end
end


function ProblemDisc:conductivity(condID)
    local conductivity = nil
    local values = self.modelMap[condID]

    if type(values) == "userdata" then
        conductivity = RichardsConductivity(self.modelMap[condID])
    end
    return conductivity
end


function ProblemDisc:saturation(satID)
    local saturation = nil
    local values = self.modelMap[satID]

    if type(values) == "userdata" then
        saturation = RichardsSaturation(self.modelMap[satID])
    end

    return saturation
end
