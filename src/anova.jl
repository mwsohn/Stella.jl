#--------------------------------------------------------------------------
# one-way or two-way analysis of variance
#--------------------------------------------------------------------------
struct ANOVAReturn
    title::String
	colnms::Vector
	array::Vector
	F::Vector
	p::Vector
	bartlett::Float64 # Bartlett's test statistic (chisq with k-1 degrees of freedom)
	bartlett_dof::Int
	bartlett_p::Float64
end

function anova(df::DataFrame,cvar::Symbol,args::Symbol...)
	if length(args) == 1
		return oneway(df,cvar,args[1])
	elseif length(args) == 2
		return towway(df,cvar,args[1],args[2])
	else
		error("More than two categorical variables are not allowed")
	end
end

#--------------------------------------------------------------------------
# oneway - return an one-way ANOVA result in ANOVAReturn
#--------------------------------------------------------------------------
function oneway(df::DataFrame,dep::Symbol,cat::Symbol)

	# establish data
	ba = completecases(df[[dep,cat]])
	dep = df[ba,dep]
	cat = df[ba,cat]

	if length(dep) == 0
		error(string(df)," does not contain usable records")
	end

	# title
	title = "One-Way Analysis of Variance"

	# N
	N = length(dep)

	# grand mean
	grandmean = mean(dep)

	# total sum of squares
	sstotal = mapreduce(x->(x-grandmean)^2,+,dep)

	# levels and number of groups
	lev = sort(collect(Set(cat)))
	k = length(lev)

	# group means and between group sum of squares
	groupmean = zeros(Float64,length(lev))
	groupvar = zeros(Float64,length(lev))
	n = zeros(Int64,length(lev))
	ssbetween = 0.
	sswithin = 0.
	for i in 1:k
		depgroup = dep[cat .== lev[i]]
		groupmean[i] = mean(depgroup)
		groupvar[i] = std(depgroup)^2
		n[i] = length(depgroup)
		ssbetween += n[i]*(groupmean[i] - grandmean)^2
		sswithin += (n[i]-1)*groupvar[i]
	end

	dfbetween = k - 1
	dfwithin = N - k
	dftotal = N - 1
	mswithin = sswithin / dfwithin
	msbetween = ssbetween / dfbetween
	mstotal = sstotal / dftotal
	F = msbetween / mswithin

	# Bartlett's chisquare of equal variance
	pooledvar = sum([(n[i] - 1)*groupvar[i] for i in 1:k])/dfwithin
	numer = dfwithin*log(pooledvar) - sum([(n[i] - 1)*log(groupvar[i]) for i in 1:k])
	denom = 1 + (sum([1/(n[i]-1) for i in 1:k]) - 1/dfwithin)/(3*dfbetween)
	bartlett = numer/denom

	return ANOVAReturn(
		title,
		["Source","SS","dof","MS","F","P-Value"],
		[
			["Between groups","Within groups","Total"],
			[ssbetween,sswithin,sstotal],
			[dfbetween,dfwithin,dftotal],
			[msbetween,mswithin,mstotal]
		],
		[F],
		[HypothesisTests.ccdf(Distributions.FDist(dfbetween,dfwithin),F)],
		bartlett,
		k - 1,
		HypothesisTests.ccdf(Distributions.Chisq(k-1),bartlett))
end

#--------------------------------------------------------------------------
# twoway - return a two-way ANOVA result in ANOVAReturn
#--------------------------------------------------------------------------
function crossmap(cat1::AbstractVector,cat2::AbstractVector)

	n = length(cat1)
	if n != length(cat2)
		error("Lengths of both vectors must be identical")
	end

	lev1 = sort(collect(Set(cat1)))
	n1 = length(lev1)
	lev2 = sort(collect(Set(cat2)))
	n2 = length(lev2)

	k=1
	d = Dict()
	for i in lev1
		for j in lev2
			d[i,j] = k
			k += 1
		end
	end

	return [d[cat1[i],cat2[i]] for i = 1:n]
end

function twoway(df::DataFrame,depvar::Symbol,cat1::Symbol,cat2::Symbol)

	# establish data
	ba = completecases(df[[depvar,cat1,cat2]])
	dep = Vector(df[ba,depvar])
	catj = Vector(df[ba,cat1])
	catk = Vector(df[ba,cat2])

	if length(dep) == 0
		error(string(df)," does not contain usable records")
	end

	# title
	title = "Two-Way Analysis of Variance"

	# N
	N = length(dep)

	# grand mean
	ybar = mean(dep)

	# levels and number of groups
	# J
	levj = sort(collect(Set(catj)))
	j = length(levj)
	ybarj = [mean(dep[catj .== i]) for i in levj]
	nj = [sum(catj .== i) for i in levj]

	# K
	levk = sort(collect(Set(catk)))
	k = length(levk)
	ybark = [mean(dep[catk .== i]) for i in levk]
	nk = [sum(catk .== i) for i in levk]

	# J x K groups
	catjk = crossmap(catj,catk)
	levjk = sort(collect(Set(catjk)))
	jk = length(levjk)
	ybarjk = [mean(dep[catjk .== i]) for i in levjk]

	# sums of squares
	ss = Vector(Float64,6)
	ss[2] = sum([nj[i]*(ybarj[i] - ybar)^2 for i in 1:j])
	ss[3] = sum([nk[i]*(ybark[i] - ybar)^2 for i in 1:k])
	ss[5] = sum([(dep[i] - ybarjk[catjk[i]])^2 for i in 1:N])
	ss[6] = sum([(dep[i] - ybar)^2 for i in 1:N])
	ss[1] = ss[6] - ss[5]
	ss[4] = ss[1] - ss[2] - ss[3]

	# degrees of freedom
	df = Vector{Int64}(6)
	df[2] = j - 1 # factor 1
	df[3] = k - 1 # factor 2
	df[4] = (j - 1)*(k - 1)
	df[1] = df[2] + df[3] + df[4]
	df[5] = N - jk
	df[6] = N - 1

	# mean squared error
	ms = ss ./ df

	# F statistics
	F = [ms[i]./ms[5] for i in 1:4]

	# P-values
	pval = [HypothesisTests.ccdf(Distributions.FDist(df[i],df[5]),F[i]) for i = 1:4]

	return ANOVAReturn(
		title,
		["Source","Partial SS","df","MS","F","P-Value"],
		[
			["Model",string(cat1),string(cat2),string(cat1,"*",cat2),"Residual","Total"],
			ss,
			df,
			ms,
		],
		F,
		pval,
		0.,
		0,
		0.)
end
#
# function anova(mout::RegressionModel; type::String = "")
# 	if in(type,["II","III"]) == false
# 		error("`type` must be either `II` or `III`")
# 	end
#
# 	if dof_residual(mout)
# 		error("residual df = 0")
# 	end
#
# 	if deviance(mout) < sqrt(eps(Float64))
# 		error("residual sum of squares is 0 (within rounding error)")
# 	end
#
# 	if typeof(mout) <: LinearModel
# 		if type == "II"
# 			return _anova_lm_II(mout,type=type)
# 		elseif type == "III"
# 			return _anova_lm_III(mout,type=type)
# 		end
# 	end
# end
#
# function _anova_lm_3(mout; type::String = "")
#
# 	# intercept
# 	intercept = mout.mf.terms.intercept
#
# 	# coef and vcov
# 	V = vcov(mout)
# 	b = coef(mout)
#
# 	# identity matrix
# 	blen = length(b)
# 	Ip = eye(Int8,blen)
#
# 	# column names
# 	Source = string.(mout.mf.terms.terms)
#
# 	# number of terms
# 	nterms = length(Source)
#
# 	# column index assignment to variables
# 	assign = mout.mm.assign
#
#
# 	# missing coefficients
# 	notaliased = !isna(coef(mout)) # find out the value returned for missing estimates
#
# 	for term in 1:nterms
# 		subs = find(x->x == term,assign)
#
# 		# hypothesis matrix
# 		hmat = Ip[subs,notaliased]
#
# 		# value
# 		hval = hmat * b
# 		hvcov = hmat * b * hmat'
#
#
#
#
#
#
# 	end
#
# end


import Base.show
Base.show(io::IO,::MIME"text/plain",a::ANOVAReturn) = show(io,a)
function show(io::IO,a::ANOVAReturn)
	width = length.(a.colnms)
	width[1] = max(width[1],maximum(length.(a.array[1])))
	width[2] = max(7,width[2],maximum(length.(Stella.strpr.(a.array[2],4)))) # 7 is minimum width
	width[3] = max(5,width[3],maximum(length.(string.(a.array[3])))) # 5 is minimum
	width[4] = max(7,maximum(length.(Stella.strpr.(a.array[4],4)))) # 7 is minimum width
	width[5] = 7 # F
	width[6] = 7 # P-value

	# title
	println(io,a.title)

	# header row
    print(io,"|  ",lpad(a.colnms[1],width[1]))
    for i=2:length(a.colnms)
        print(io,"  |  ",lpad(a.colnms[i],width[i]))
    end
    print(io,"  |\n")

	# line
	Stella.hline(io,width.+2)

	# array output
	len = length(a.array[1])
	for i = 1:len-1
		print(io,"|  ",lpad(a.array[1][i],width[1]))
		for j=2:4
			print(io,"  |  ",lpad(j == 3 ? a.array[j][i] : strpr(a.array[j][i],4),width[j]))
		end
		if i < len-1
			print(io,"  |  ",lpad(strpr(a.F[i],5),width[5]))
			print(io,"  |  ",lpad(strpr(a.p[i],5),width[6]))
			print(io,"  |")
		else
			print(io,"  |  ",repeat(" ",width[5]),"  |  ",repeat(" ",width[6]),"  |")
		end
		print(io,"\n")
	end

	# line
	Stella.hline(io,width.+2)

	# total
	print(io,"|  ",Stella.lpad(a.array[1][len],width[1]))
	for j=2:4
		print(io,"  |  ",Stella.lpad(j == 3 ? a.array[j][len] : strpr(a.array[j][len],4),width[j]))
	end
	print(io,"  |  ",repeat(" ",width[5]),"  |  ",repeat(" ",width[6]),"  |")
	print(io,"\n\n")
end
