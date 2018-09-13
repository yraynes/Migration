#!/contrib/projects/julia/julia
using Distributions

type Lineage
  fitness::Float64
  size::Int64
  state::Vector{Float64}
  mut_rate::Float64
end

function average_fitness(population::Dict{Array{Float64,1}, Lineage})
#calculate average fitness of the population. also outputs population size
  popN = Int64
  popN = 0
  popW = Float64
  popW = 0
  for line in values(population)
    popN+=line.size
    popW+=line.size * line.fitness
  end
  popW = popW/popN
  return popW, popN
end

function average_meta_fitness(metapopulation::Array{Dict{Array{Float64,1}, Lineage}})
    popN = Int64
    popN = 0
    popW = Float64
    popW = 0
    for p = 1:length(metapopulation)
        for line in values(metapopulation[p])
          popN+=line.size
          popW+=line.size * line.fitness
        end
    end
    popW = popW/popN
    return popW, popN
end

function pop_size(population::Dict{Array{Float64,1}, Lineage})
  popN = Int64
  popN = 0
  for line in values(population)
    popN+=line.size
  end
  return popN
end

function assay_mutation_rate(population::Dict{Array{Float64,1}, Lineage})
#calculates mutator frequency
  popN = Int64
  popN = 0
  mutN = Int64
  mutN = 0
  for line in values(population)
    popN+=line.size
    if line.mut_rate>1 mutN+=line.size end
  end
  mutF = Float64
  mutF = mutN/popN
  return mutF
end

function assay_meta_mutation_rate(metapopulation::Array{Dict{Array{Float64,1}, Lineage}})
#calculates mutator frequency
  popN = Int64
  popN = 0
  mutN = Int64
  mutN = 0
  mutF = Float64
  for p = 1:length(metapopulation)
      for line in values(metapopulation[p])
        popN+=line.size
        if line.mut_rate>1 mutN+=line.size end
      end
  end
  mutF = mutN/popN
  return mutF
end


function fitness_calc(state::Array{Float64,1})
#calculates fitness from the lineage state = sum of all mutational effects - the mutator strength (100) - oriT marker (-100)
  w = Float64
  w = 1 + sum(state) - state[1]
  return w
end

function avail_sites(state::Array{Float64,1})
  sites = Int64[]
  for i = 1:length(state)
    if state[i] == 0 push!(sites,i) end
  end
  return sites
end

function mutate_population(population::Dict{Array{Float64,1}, Lineage}, Ub::Float64, sb::Float64, Ud::Float64, sd::Float64, mut_multiplier::Float64, mode)
  new_population = Dict{Array{Float64,1}, Lineage}()
  #initialize new population dictionary
  for (key, line) in population
	if mode == 1
    mutations = rand(Poisson(line.mut_rate * (Ud + Ub) * line.size))
	elseif mode == 2
	mutations = 0
	end
    if mutations > line.size mutations = line.size end
    bmutations = rand(Binomial(mutations, Ub/(Ub+Ud)))
    dmutations = mutations - bmutations

    open_sites =   avail_sites(line.state)
    #generates a list of indexes of all loci that are available for mutation
    mutation_positions = open_sites[rand(1:end,(dmutations+bmutations))] #generates positions of new mutations by randomly picking indexes from open_sites

      #for each new mutation, copies the state, adds mutation and generates the key.
      #check if it is alrady in the new_population. if it is, add 1 individual to it, if it's not - make another lineage with
      #the key and size 1
    c=1
    for i = 1:bmutations
      #assign beneficial mutations
      newstate = copy(line.state)#
      newstate[mutation_positions[c]] += sb
      c+=1
      if newstate in keys(new_population) new_population[newstate].size+=1
      else
        newfit = 1 + sum(newstate[2:end])
        new_population[newstate] = Lineage(newfit, 1, newstate, newstate[1])
      end
    end

    for i = 1:dmutations
      #assign deleterious mutations
      newstate = copy(line.state)#
      newstate[mutation_positions[c]] += -sd
      c+=1
      if newstate in keys(new_population) new_population[newstate].size+=1
      else
        #newfit = Float64
        #newfit::Float64 = 1 + sum(newstate)-newstate[1]
        newfit::Float64 = 1 + sum(newstate[2:end])
        new_population[newstate] = Lineage(newfit, 1, newstate, newstate[1])
      end
    end
    new_size = line.size - bmutations - dmutations
    if  new_size > 0
      if key in keys(new_population) new_population[key].size+=new_size
      else new_population[key] = Lineage(line.fitness, new_size, line.state, line.mut_rate)
      end
    end
  end
  return new_population
end
#end

# import Distributions.isprobvec
# isprobvec(p::Vector{Float64}) = true

function wright_fisher_reproduction(population::Dict{Array{Float64,1}, Lineage}, N0::Int64, Nnew::Int64, popw::Float64)
  proby_list = Float64[] #array of probabilities
  lineage_list = Lineage[] #initialize lineage liste
  for (key,line) in population
    push!(proby_list, line.size/N0 * line.fitness/popw)
    #representation of a lineage in the next generations depends on its frequency and relative fitness
    push!(lineage_list, line)
  end
  new_counts_list = rand(Multinomial(Nnew, proby_list))
  new_population = Dict{Array{Float64,1}, Lineage}()
  for i = 1:length(new_counts_list)
    if new_counts_list[i] > 0
      new_population[lineage_list[i].state] = Lineage(lineage_list[i].fitness, new_counts_list[i], lineage_list[i].state, lineage_list[i].mut_rate)
    end
  end
  return new_population
end

function combine_pops(population1::Dict{Array{Float64,1}, Lineage},population2::Dict{Array{Float64,1}, Lineage})
  new_population = Dict{Array{Float64,1}, Lineage}()
  for (key,line) in population1
    if key in keys(new_population) new_population[key].size+=line.size
    else new_population[key] = Lineage(line.fitness, line.size, line.state, line.mut_rate)
    end
  end
  for (key,line) in population2
    if key in keys(new_population) new_population[key].size+=line.size
    else new_population[key] = Lineage(line.fitness, line.size, line.state, line.mut_rate)
    end
  end
  return new_population
end

function bottleneck(population::Dict{Array{Float64,1}, Lineage}, N0::Int64, Nnew::Int64)
  #not for the HGT project. subjects populations to random bottlenecks
  proby_list = Float64[] #array of probabilities
  lineage_list = Lineage[]
  for (key,line) in population
    push!(proby_list, line.size/N0)
    push!(lineage_list, line)
  end
  new_counts_list = rand(Multinomial(Nnew, proby_list))
  new_population = Dict{Array{Float64,1}, Lineage}()
  for i = 1:length(new_counts_list)
    if new_counts_list[i] > 0
      new_population[lineage_list[i].state] = Lineage(lineage_list[i].fitness, new_counts_list[i], lineage_list[i].state, lineage_list[i].mut_rate)
    end
  end
  return new_population
end

function bottleneck2(population::Dict{Array{Float64,1}, Lineage}, N0::Int64, Nnew::Int64)
  #random bottleneck without replacement, returns remainder of population and the sample called new_population
  key_list = Array{Array{Float64,1}}(N0)
  i::Int64=1
  for (key,line) in population
    for j = 1:line.size
      key_list[i] = key
      i+=1
    end
  end
  key_list = sample(key_list,Nnew,replace=false)
  new_population = Dict{Array{Float64,1}, Lineage}()
  for newstate in key_list
      if newstate in keys(new_population) new_population[newstate].size+=1
      else
        newfit = 1 + sum(newstate[2:end])
        new_population[newstate] = Lineage(newfit, 1, newstate, newstate[1])
      end
      population[newstate].size -= 1
  end
  return population, new_population
end

function one_generation(population::Dict{Array{Float64,1}, Lineage}, mode::Int64, growth::Int64, mutator_strength::Float64, Ub::Float64, sb::Float64, Ud::Float64, sd::Float64)
    popw, pop_N = average_fitness(population)
    population = wright_fisher_reproduction(population, pop_N, pop_N*growth, popw)
    population = mutate_population(population, Ub, sb, Ud, sd, mutator_strength, mode)
    mut_f = assay_mutation_rate(population)
    popw, pop_N = average_fitness(population)
return population, popw, pop_N, mut_f
end

function one_meta_generation(metapopulation::Array{Dict{Array{Float64,1}, Lineage}}, mode::Int64, growth::Int64, mutator_strength::Float64, Ub::Float64, sb::Float64, Ud::Float64, sd::Float64)
    for i = 1:length(metapopulation)
        popw, pop_N = average_fitness(metapopulation[i])
        metapopulation[i] = wright_fisher_reproduction(metapopulation[i], pop_N, pop_N*growth, popw)
        metapopulation[i] = mutate_population(metapopulation[i], Ub, sb, Ud, sd, mutator_strength, mode)
    end
return metapopulation
end

function combine_metapop(metapopulation::Array{Dict{Array{Float64,1}, Lineage}})
  new_population = Dict{Array{Float64,1}, Lineage}()
  for i = 1:length(metapopulation)
      for (key,line) in metapopulation[i]
        if key in keys(new_population)
            new_population[key].size+=line.size
        else
            new_population[key] = Lineage(line.fitness, line.size, line.state, line.mut_rate)
        end
      end
  end
  return new_population
end

function meta_migration_1(metapopulation::Array{Dict{Array{Float64,1}, Lineage}}, pop_Ni::Int64)
    #at each migration event all demes are combined and re-drawn from the mix without replacement
    mix = combine_metapop(metapopulation)
    new_metapopulation = Array{Dict{Array{Float64,1}, Lineage}}(length(metapopulation))
    for i = 1:length(metapopulation)
        mix, new_metapopulation[i] = bottleneck2(mix, pop_size(mix), pop_Ni)
    end
return new_metapopulation
end

function meta_migration_2(metapopulation::Array{Dict{Array{Float64,1}, Lineage}}, pop_Ni::Int64, P::Float64)
    mix = Dict{Array{Float64,1}, Lineage}()
    migrants_out = rand(Poisson(pop_Ni*P),length(metapopulation))
    for i = 1:length(metapopulation)
        metapopulation[i], add_to_mix  = bottleneck2(metapopulation[i], pop_size(metapopulation[i]), min(pop_Ni,migrants_out[i]))
        mix = combine_pops(mix, add_to_mix)
    end
    for i = 1:length(metapopulation)
        migrants_in = pop_Ni-pop_size(metapopulation[i])
        mix, add_to_pop = bottleneck2(mix, pop_size(mix), migrants_in)
        metapopulation[i] = combine_pops(metapopulation[i], add_to_pop)
    end

return metapopulation
end

function simulate(pop_Ni, T, neutral,Pc, B, r)
  #parameters defined first
  const sb = 0.1
  const sd = 0.1
  const Ub = 0.000001
  const Ud = 0.0001
  const mutator_strength = 100.0

  metapopulation = Array{Dict{Array{Float64,1}, Lineage}}(Pc)
  wt_state = zeros(Float64,100)
  wt_state[1] = 1.0
  mut_state = deepcopy(wt_state)
  mut_state[1] = mutator_strength
  mut_state[2] = B
  for i = 1:Pc
      init_mut_N = 100
      init_wt_N = pop_Ni-init_mut_N
      metapopulation[i] = Dict{Array{Float64,1}, Lineage}()
      metapopulation[i][wt_state]  = Lineage(1,init_wt_N, wt_state, 1.0)
      metapopulation[i][mut_state] = Lineage(1+B,init_mut_N, mut_state, mutator_strength)
  end


  mut_av = assay_meta_mutation_rate(metapopulation)


  generations = Int64
  generations = 1
  mw, mN = average_meta_fitness(metapopulation)
  #println("At the beginning: ", generations, " Mut frequency: ", mut_av, " Meta Fitness: ", mw, " Population size: ", mN)
  growth::Int64=1
  #initialize generations counter
  while 0.0<mut_av<1.0

    metapopulation = one_meta_generation(metapopulation,  neutral, growth, mutator_strength, Ub, sb, Ud, sd)

    if Pc > 1
        if generations/T == Int(round(generations/T))
            metapopulation = meta_migration_2(metapopulation, pop_Ni, r)
        end
        mut_av = assay_meta_mutation_rate(metapopulation)
    elseif Pc == 1
        mut_av = assay_mutation_rate(metapopulation[1])
    end

    generations+=1
  end
  return mut_av ==1.0, generations, average_meta_fitness(metapopulation)[1]
end

job_id = ARGS[1]

neutral = 1#switch to make the mutator just a neutral allele (used to test simulation)
justB = 0#make a mutator inherently beneficial (used to test simulation)

Pc = 6
n = 200

for r = [0.01, 0.1, 0.2, 0.5, 1.0] #percentage of each deme donated to common pool during migration
    for T = [1, 10, 20, 50, 75, 100, 150, 200, 250] #time between migration events
          time_to_mut= Float64[]
          #time_to_mut list records all times of successful mutator hitchhiking
          for run = 1:2000
            output = simulate(n, T, neutral, Pc, justB, r)
            if output[1]
              push!(time_to_mut, output[2])
            else
              push!(time_to_mut, 0)
            end
            end
          outfile = open(string(job_id,"time_to_mut.csv"), "a")
          write(outfile, join(time_to_mut, ","), "\n")
          close(outfile)
    end
end
