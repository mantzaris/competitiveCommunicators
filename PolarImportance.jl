using Plots
using StatsBase


function mainAlex(NN=20,TT=8000,pb=0.01,cb=1,cr=4)
    Ratios = Array(Float64,8,NN)
    Importances = Array(Float64,8,NN)
    filename = "ScatterPlots/x3ImportTest_";
    AA = dynamicCommunicators(NN,TT,pb,cb,cr,filename,Ratios,Importances)
    degOut = totalOut(AA,NN)
    Cbroad = dynamicCentrality(AA,NN,TT)
    scatterPlot(AA,size(AA)[1],degOut,Cbroad)
    savefig("ScatterPlots/x3ImportTest_Final");
    TimeSeries(Ratios,1);
    savefig("ScatterPlots/Cbroad_Evolving")
    TimeSeries(Importances, 0);
    savefig("ScatterPlots/Importance_Evolving")
    rankTimeSeries(Ratios,NN)
    

   # AA2 = ringNet()
   # degOut = totalOut(AA2,21)
   # Cbroad = dynamicCentrality(AA2,21,3)
   # scatterPlot(AA2,size(AA2)[1],degOut,Cbroad)
end

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
    title!(string("Number of Nodes = 20  -  Importance X3  -  Network Evolving - Polarized")) # UPDATE THIS BEFORE YOU PLOT AGAIN
end

function scatterPlot(AA,NN,degOut,Cbroad)
    labels = [string(x) for x in 1:NN]
    scatter(degOut,Cbroad,markersize=1,series_annotations =labels,grid=false)
    xlabel!("total deg out")
    ylabel!("broadcast centrality")
    title!(string("node num=",NN))
    
end

function TimeSeries(Ratios, n)
    plot(Ratios, label = [string("Node ", x) for x in 1:20])
    if(n == 1)
    	xlabel!("Time")
    	ylabel!("Responses Triggered")
    else
	xlabel!("Time")
    	ylabel!("Importance")
    end
    	title!(string("Number of Nodes = 20"))

end


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


function dynamicCommunicators(NN,TT,pb,cb,cr,filename,Ratios,Importances)
    imp = [x^3 for x in 1:NN]
    initial_imp = imp;
    midpoint_imp = ones(NN,1)
    i = 1;
    println(imp)
    println("that was imp vector now normal")
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

	if(tt == 8000)
	   Mask1 = [imp[x]*((-1)^x) for x in 1:NN]
   	   Mask2 = [imp[x]*((-1)^(x+1)) for x in 1:NN] 
	end

	if(tt == (TT/2))
	    midpoint_imp = imp
	end
	if(tt%(TT/8) == 0) # Print Scatter plot 4 times evenly spaced
	    degOut = totalOut(AA,NN);
    	    Cbroad = dynamicCentrality(AA,NN,tt);
            scatterPlot(AA,size(AA)[1],degOut,Cbroad);
            savefig(string(filename,(string(i))));
	    Ratios[i,:] = (.2*(imp-initial_imp));
            Importances[i,:] = (imp-initial_imp)
	    i+=1;
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

                if(ii%2 == 0)
                    totalImportance = sum([msgsTo_ii[x]*Mask1[x] for x in 1:NN]); # Even Node gets Odd Importance
		end
 		else
		    totalImportance = sum([msgsTo_ii[x]*Mask2[x] for x in 1:NN]); # Odd Node gets Even Importance
                end

                #println(msgsTo_ii)
                #println(imp)
                #println(sum([msgsTo_ii[x]*imp[x] for x in 1:NN]))
                
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
                end
            end
        end
                
    end
    for ii in 1:NN
	println(ii, " Initial Importance[ii] = ", initial_imp[ii], " Final Importance[ii] = ", imp[ii], " Growth = ", imp[ii]-initial_imp[ii], 		" Responses Generated = ", (.2)*(imp[ii]-initial_imp[ii]))
    end  
    println("Spearman Correlation Initial -> Mid = ", corspearman(initial_imp, midpoint_imp), " initial -> final = ", corspearman(initial_imp, 	   imp)) 
    println("Kendall Correlation Initial -> Mid = ", corspearman(initial_imp, midpoint_imp), " initial -> final = ", corspearman(initial_imp, imp))     
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



