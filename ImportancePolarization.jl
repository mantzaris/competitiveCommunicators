using Plots
using StatsBase
using PlotRecipes
gr()

# NOTES FOR LATER
##################################################################################################

# 1) Keep the cbroad/degout graph.
# 2) Find a way to show difference cbroad/averagecbroad or maybe divide by magnitude of total degout.

##################################################################################################

function mainAlex(NN=20,TT=16000,pb=0.01,cb=1,cr=4)

    #Passing all these in from outside because plotting doesn't work otherwise
    Responses = Array(Float64,8,NN)    # Total Responses from each Node at 8 Time Points
    Ratios = Array(Float64,8,NN)       # Cbroad/degOut for each node at 8 Time Points
    Importances = Array(Float64,8,NN)  # Importance Growth Final - Initial for each node
  #  MatSum = Array(Float64,NN,NN)      # Aggregates messages from nodes by += every time point
    MatSumMid = Array(Float64,NN,NN,4) # Sum of all messages sent from all nodes at 4 time points 1,2 Pre Polar - 3,4 Post Polar
    BB = Array(Float64,NN,NN,2)        # dynamicCentrality Matrix at time = 8000 and time = 16000
    
    # Main Procedure
    AA = dynamicCommunicators(NN,TT,pb,cb,cr,Ratios,Importances,Responses,MatSumMid,BB);

    # Create Plots # Figure out how to Plot correctly
    plot(BB[:,:,1])#label=[string("Node ", x) for x in 1:20],grid=false,background_color=RGB(.2,.2,.2))



    
end

# Scatter Plot Function
function scatterPlot(AA,NN,degOut,Cbroad)
    labels = [string(x) for x in 1:NN]
    scatter(degOut,Cbroad,markersize=1,series_annotations=labels,grid=false)
    xlabel!("Total Deg Out")
    ylabel!("Broadcast Centrality")
    title!(string("Node Num = ",NN))
    
end

# Function to plot Time Series of Rank
function rankTimeSeries(Ratios,NN)
    Rank = Array(Float64,8,NN)

    for ii in 1:8  # Rank[ii,jj] = Num nodes lower than node ii,jj at timepoint (ii*1000)
        for jj in 1:NN
	    Rank[ii,jj] = sum([Ratios[ii,jj] < x for x in Ratios[ii,:]])
	end
    end

    plot(Rank,label=[string("Node ", x) for x in 1:20],grid=false,background_color=RGB(.2,.2,.2))
    xlabel!("Time")
    ylabel!("Rank")
    title!(string("Number of Nodes = 20  -  Importance X3  -  Network Evolving")) # UPDATE THIS BEFORE YOU PLOT AGAIN
end

# Function to Plot Time Series of Ratio - Cbroad/degOut
function TimeSeries(Ratios, n)

    plot(Ratios,label=[string("Node ", x) for x in 1:20],grid=false,background_color=RGB(.2,.2,.2))
    # NEEDS LABELS AND TITLE DONT FORGET BEFORE YOU PLOT

end

# Produces Matrix based on Set Data
function ringNet()
    AA2 = zeros(21,21,3)
    AA2[1,7,1] =1;AA2[7,1,1] =1;AA2[2,7,1] =1;AA2[7,2,1] =1;
    AA2[21,20,1] =1;AA2[20,21,1] =1;AA2[21,19,1] =1;AA2[19,21,1] =1;
    AA2[2,8,2]=1;AA2[8,2,2]=1;AA2[3,8,2]=1;AA2[8,3,2]=1;
    AA2[20,18,2]=1;AA2[18,20,2]=1;AA2[20,17,2]=1;AA2[17,20,2]=1;
    AA2[19,16,2]=1;AA2[16,19,2]=1;
    AA2[4,10,3]=1;AA2[10,4,3]=1;AA2[5,10,3]=1;AA2[10,5,3]=1;
    AA2[18,15,3]=1;AA2[15,18,3]=1;AA2[18,14,3]=1;AA2[14,18,3]=1;
    AA2[17,14,3]=1;AA2[14,17,3]=1;AA2[17,13,3]=1;AA2[13,17,3]=1;
    AA2[16,13,3]=1;AA2[13,16,3]=1;AA2[16,12,3]=1;AA2[12,16,3]=1;
    return AA2
end


function dynamicCommunicators(NN,TT,pb,cb,cr,Ratios,Importances,Responses,MatSumMid,BB)

    imp = [x^3 for x in 1:NN]
    initial_imp = imp;
    midpoint_imp = ones(NN,1)
    response = zeros(NN,1)
    i = 1;


    AA = zeros(NN,NN,TT)
    #INIT DATA>>
    for ii in 1:NN
        if(rand(1)[1] <= pb)
            tmp = deleteat!(collect(1:NN),ii)
            destinationN = tmp[rand(1:end)]            
            AA[ii,destinationN,1] = 1
        end

    end
    #DATA LOOP: 1st put basal edges into <present> / 2nd put response from past to now
    for tt in 2:(TT)
	
	# Capture MatSum from Mid of non Polarized Run
	if(tt == (TT/4))
	    # midpoint_imp = imp # Legacy from calculating Tau/Rho
  	    MatSumMid[:,:,1] = sum(AA,3);
	end
 
	# Capture BB & MatSum from end of non Polarized Run
	if(tt == (TT/2))
	    BB[:,:,1] = dynamicCentrality2(AA,NN,TT);
	    MatSumMid[:,:,2] = sum(AA,3);
	end

	# Capture MatSum from Mid of Polarized Run
	if(tt == ((3*TT)/4))
	    MatSumMid[:,:,3] = sum(AA,3);
	end

	# Capture BB & MatSum from end of Polarized Run
	if(tt == TT)
	    BB[:,:,2] = dynamicCentrality2(AA,NN,TT);
	    MatSumMid[:,:,4] = sum(AA,3);
	end

	# Legacy from Plotting Cbroad/degOut, Importances, & Responses
#	if(tt%(TT/4) == 0) # Calculate Responses, Importances, Ratio Cbroad/degOut at 8 evenly spaced time points
#	    degOut = totalOut(AA,NN);
#    	    Cbroad = dynamicCentrality(AA,NN,tt);
#	    Responses[i,:] = response;
#           Importances[i,:] = (imp-initial_imp)
#	    Ratios[i,:] = Cbroad./degOut;
#	    i+=1;
#       end

	# After 8000 timepoints Polarization is introduced 
        if(tt == 8001)
	   imp = [imp[x]*((-1)^x) for x in 1:NN]
	end 

        #BASAL LOOP
        for ii in 1:NN
            if(rand(1)[1] <= pb)
                tmp = deleteat!(collect(1:NN),ii)
                destinationN = tmp[rand(1:end)]            
                AA[ii,destinationN,tt] = 1#FIRE NOW
            end

        end
        #RESPONSE LOOP
        for ii in 1:NN
            msgsTo_ii = AA[:,ii,tt-1]#RESPOND TO PREVIOUS STATE (tt-1)          
            if(sum(AA[:,ii,tt-1]) >= 1)
                    
                totalImportance = sum([msgsTo_ii[x]*imp[x] for x in 1:NN]);                
                r_nNumerator = totalImportance
                r_nDenominator = 1 + (findmax(imp)[1] * sum(AA[:,ii,tt-1]))
                r_next = r_nNumerator / r_nDenominator
                println(r_next)
                if(r_next >=  rand(1)[1])#if so generate cr links

                    tmp = deleteat!(collect(1:NN),ii)
                    destinationNodes = tmp[randperm(length(tmp))[1:cr]]
                    AA[ii,destinationNodes,tt] = 1#FIRE NOW
                    println("fire!")
		    imp = ((5)*msgsTo_ii + imp);
  		    response = response + msgsTo_ii;
                end
            end
        end
                
    end
    return AA
end


function totalOut(AA,NN)
    coltmp = zeros(NN,1)
    outDegVec = sum(sum(AA,3),2)
    coltmp[:] = outDegVec[:]
    return coltmp
end


function dynamicCentrality(AA,NN,TT)
    println(":")
    daynorm = zeros(TT,1);
    for ii in 1:TT
        Atemp = AA[:,:,ii]
        daynorm[ii] = findmax(abs(eigvals(Atemp)))[1]
    end    
    alpha = 0.9*(1/findmax(daynorm)[1]);#90percent of the max
    println(alpha)
    #alpha = 0.15;
    BB = eye(NN,NN);
    for kk = TT:-1:1
        Aday = AA[:,:,kk]
        #println(Aday)
        BBtmp = inv(eye(NN,NN) - alpha*Aday)#\BB
        #BBtmp = (eye(NN,NN) - alpha*Aday) \ BB
        BBtmp = abs(BBtmp)
        #BB = BBtmp / norm(BBtmp)
        BB = (BB*BBtmp) / norm(BB)
    end
    #sum across columns
    Cbroad = sum(BB',2)    
    return Cbroad
end

function dynamicCentrality2(AA,NN,TT)
    println(":")
    daynorm = zeros(TT,1);
    for ii in 1:TT
        Atemp = AA[:,:,ii]
        daynorm[ii] = findmax(abs(eigvals(Atemp)))[1]
    end    
    alpha = 0.9*(1/findmax(daynorm)[1]);#90percent of the max
    println(alpha)
    #alpha = 0.15;
    BB = eye(NN,NN);
    for kk = TT:-1:1
        Aday = AA[:,:,kk]
        #println(Aday)
        BBtmp = inv(eye(NN,NN) - alpha*Aday)#\BB
        #BBtmp = (eye(NN,NN) - alpha*Aday) \ BB
        BBtmp = abs(BBtmp)
        #BB = BBtmp / norm(BBtmp)
        BB = (BB*BBtmp) / norm(BB)
    end
    return BB;   
end


