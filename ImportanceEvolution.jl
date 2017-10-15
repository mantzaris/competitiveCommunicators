using Plots
using StatsBase
gr()

function main(NN=40,TT=3650,pb=0.01,cb=1,cr=4)
    	
     imp1 = ones(NN,1)*10; # Constant

    # Extra Importantces for additional runs
 #   imp2 = imp1;  
 #   imp2[NN] += 90; #Constant with bigger final
 #   imp3 = [x for x in 1:NN]  # Linear
     imp4 = [e^x for x in 1:NN] # exponential
 #   imp5 = rand(1:100, NN)
 #   imp5[NN] = 1000
 #   imp5[20] = 3150
 #   imp5[35] = 3150
 #   imp5[27] = 4200

    names = ["ScatterPlots/All_10s_Imp_", "ScatterPlots/10s_with_one_100_", "ScatterPlots/Linear_Imp_", "ScatterPlots/Exponential_Imp_", "ScatterPlots/Random_Imp_1-100_One_1000_"]

    # Exponential Importance Run
    AA = dynamicCommunicators(NN,TT,pb,cb,cr,imp4,names[4])
    degOut = totalOut(AA,NN)
    Cbroad = dynamicCentrality(AA,NN,TT)
    scatterPlot(AA,size(AA)[1],degOut,Cbroad)
    savefig("ScatterPlots/Exponential_Imp_Final")
    

    # Extra Calls for Extra Importance Functions
  # 10s with one 100 Imp Run
  #  AA = dynamicCommunicators(NN,TT,pb,cb,cr,imp2,names[2])
  #  degOut = totalOut(AA,NN)
  #  Cbroad = dynamicCentrality(AA,NN,TT)
  #  scatterPlot(AA,size(AA)[1],degOut,Cbroad)
  #  savefig("ScatterPlots/10s_with_one_100_Final")
    
  # Linear Imp Run
  #  AA = dynamicCommunicators(NN,TT,pb,cb,cr,imp3,names[3])
  #  degOut = totalOut(AA,NN)
  #  Cbroad = dynamicCentrality(AA,NN,TT)
  #  scatterPlot(AA,size(AA)[1],degOut,Cbroad)
  #  savefig("ScatterPlots/Linear_Imp_Final")

  # Exponential Imp Run
  #  AA = dynamicCommunicators(NN,TT,pb,cb,cr,imp4,names[4])
  #  degOut = totalOut(AA,NN)
  #  Cbroad = dynamicCentrality(AA,NN,TT)
  #  scatterPlot(AA,size(AA)[1],degOut,Cbroad)
  #  savefig("ScatterPlots/Exponential_Imp_Final")

  # All random below 100 One 1000 Imp Run
  #  AA = dynamicCommunicators(NN,TT,pb,cb,cr,imp5,names[5])
  #  degOut = totalOut(AA,NN)
  #  Cbroad = dynamicCentrality(AA,NN,TT)
  #  scatterPlot(AA,size(AA)[1],degOut,Cbroad)
  #  savefig("ScatterPlots/Random_Imp_1-100_One_1000_Final")

   # AA2 = ringNet()
   # degOut = totalOut(AA2,21)
   # Cbroad = dynamicCentrality(AA2,21,3)
   # scatterPlot(AA2,size(AA2)[1],degOut,Cbroad)
end


function scatterPlot(AA,NN,degOut,Cbroad)
    labels = [string(x) for x in 1:NN]
    scatter(degOut,Cbroad,markersize=1,series_annotations =labels,grid=false)
    xlabel!("total deg out")
    ylabel!("broadcast centrality")
    title!(string("node num=",NN))
    
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


function dynamicCommunicators(NN=40,TT=50,pb=0.01,cb=1,cr=4, imp = [x for x in 1:NN], name = "null")
    initial_imp = imp; # Capture Initial Importance
    midpoint_imp = ones(NN,1) # Initialize Midpoint
    AA = zeros(NN,NN,TT)
    i=1;
    n = .01*findmax(imp)[1]; # Maximum importance 1% of maximum for adding on successful responses
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
	if(tt%850 == 0) # Print Scatter plot 4 times evenly spaced
	    degOut = totalOut(AA,NN);
    	    Cbroad = dynamicCentrality(AA,NN,TT);
            scatterPlot(AA,size(AA)[1],degOut,Cbroad);
            savefig(string(name,(string(i))));
	    i+=1;
        end
	if(tt == (3650/2)) # grab importance halfway through
	    println("it worked")
     	    midpoint_imp = imp;
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
            #totalImportance = sum(msgsTo_ii' .* imp)
            if(sum(AA[:,ii,tt-1]) >= 1)
            #   println(sum(msgsTo_ii .* imp  ))
            #   println((msgsTo_ii .* imp  ))                            
                totalImportance = sum( msgsTo_ii .* imp)
                r_nNumerator = totalImportance
                r_nDenominator = 1 + (findmax(imp)[1] * sum(AA[:,ii,tt-1]))
            #   println( r_nNumerator / r_nDenominator )
                r_next = r_nNumerator / r_nDenominator

                if(r_next >  rand(1)[1])#if so generate cr links

                    tmp = deleteat!(collect(1:NN),ii)
                    destinationNodes = tmp[randperm(length(tmp))[1:cr]]
                    #AA[ii,destinationNodes,tt] = 1#FIRE NOW
		    imp = (n*msgsTo_ii + imp) # Add 1% of the maximum importance to each imp of each node that triggered successful response

               end
            end
        end
                
    end
    for ii in 1:NN
	println(ii, " Initial Importance[ii] = ", initial_imp[ii], " Final Importance[ii] = ", imp[ii], " Growth = ", imp[ii]-initial_imp[ii], " Responses Generated = ", (1/n)*(imp[ii]-initial_imp[ii]))
    end  
    println("Spearman Correlation Initial -> Mid = ", corspearman(initial_imp, midpoint_imp), " initial -> final = ", corspearman(initial_imp, imp)) 
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
    #alpha = 0.9*(1/findmax(daynorm)[1]);#90percent of the max
    #println(alpha)
    alpha = 0.65;
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

