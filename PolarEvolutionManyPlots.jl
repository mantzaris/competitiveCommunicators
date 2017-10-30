#Pkg.checkout("PlotRecipes")
using StatsBase
using Plots
using PlotRecipes
using LightGraphs
using StatPlots
#using GraphPlot
gr(size=(1000,1000))


#pyplot(size=(1000,1000))

function mainAlex(NN=20,TT=24000,pb=0.01,cb=1,cr=4)
    
    MessagesFired = Array{Float64}(24,NN)
    ResponsesTriggered = zeros(NN,NN,2)
    
    AA = dynamicCommunicators(NN,TT,pb,cb,cr,MessagesFired,ResponsesTriggered)
    degOut = totalOut(AA,NN)
    Cbroad = dynamicCentrality(AA,NN,TT)
#    scatterPlot(AA,size(AA)[1],degOut,Cbroad)
#    rankTimeSeries(MessagesFired,NN);
#    TimeSeries(MessagesFired,NN);

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

function multiHistogram(data,file,title)
    labels = ["Pre-Polarization", "Post-Polarization"]
    groupedbar(data, labels=labels, nbins=20)
    xlabel!("Nodes")
    ylabel!("Broadcast Power Relative to Weakest Node")
    title!(title)
    savefig(file)
end

function graphVisualization(G,NN,file)
    names = [string(x) for x in 1:NN]
  #  g = DiGraph(G)
    graphplot(G,nodesize=22,names=names)
    savefig(file)
    println("Graph Visualization")
    
end

function heatMap(AA,file,flag)
        xs = [string(x) for x in 1:20]
	ys = [string(x) for x in 1:20]
	heatmap(xs,ys,AA);
	title!(string("Send vs Recieve ", (flag==1) ? "Non Polar" : "Polarized"))
	savefig(file)
	println("Heat Map")
end

function rankTimeSeries(MessagesFired,NN)
    Rank = Array{Float64}(16,NN)

    for ii in 1:16
    # Rank[ii,jj] = Num nodes lower than node ii,jj at timepoint (ii*1000)
        for jj in 1:NN
	    Rank[ii,jj] = sum([MessagesFired[ii,jj] < x for x in MessagesFired[ii,:]])
	end
    end
    plot(Rank,label=[string("Node ", x) for x in 1:20],grid=false,background_color=RGB(.2,.2,.2))
    xlabel!("Time")
    ylabel!("Rank")
    title!(string("Number of Nodes = 20 - Rank vs Time - Polar Importance"))
    savefig("Plots/Rank_Importance_Evolving_Polar");
end

function TimeSeries(Data,NN,file,title)
   
    plot(Data,label=["Power Top/Bot"],grid=false,background_color=RGB(.2,.2,.2))
    xlabel!("Time")
    ylabel!("Average Power Top Nodes over Bottom Nodes")
    title!(title)
    savefig(file);
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


function dynamicCommunicators(NN,TT,pb,cb,cr,MessagesFired,ResponsesTriggered)

    imp = [x^3 for x in 1:NN]
    Mask1 = [(-1)^x for x in 1:NN]
    Mask2 = [(-1)^(x+1) for x in 1:NN]
    Messages = Array{Float64}(20)
    data = Array{Float64}(20,2)
    Avg = Array{Float64}(24)
    d = 1
    f = 1
    AA = zeros(NN,NN,TT)
    i=1;

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
        #println(string("tt=",tt))

        # Collect Messages Fired to Calculate Rank
        if(tt%1000 == 0)
            MessagesFired[i,:] = Messages;
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
	    
	        if(tt <= 8000)
                    totalImportance = sum([msgsTo_ii[x]*imp[x] for x in 1:NN])
		else
		    if(ii%2==0)
                       totalImportance = sum([msgsTo_ii[x]*imp[x]*Mask1[x] for x in 1:NN])
                    else
                       totalImportance = sum([msgsTo_ii[x]*imp[x]*Mask2[x] for x in 1:NN])
	            end
		end
		
                r_nNumerator = totalImportance
                r_nDenominator = 1 + (findmax(imp)[1] * sum(AA[:,ii,tt-1]))
                r_next = r_nNumerator / r_nDenominator
              #  println(r_next)
                if(r_next >=  rand(1)[1])#if so generate cr links

                    tmp = deleteat!(collect(1:NN),ii)
                    destinationNodes = tmp[randperm(length(tmp))[1:cr]]
                    AA[ii,destinationNodes,tt] = 1#FIRE NOW
		    imp =5*(msgsTo_ii) + imp;
		    Messages = Messages + msgsTo_ii;
		    if(tt < 8000)
		        ResponsesTriggered[:,ii,1] += msgsTo_ii;	
		    else
		        ResponsesTriggered[:,ii,2] += msgsTo_ii;
		    end
                end
            end
        end
    end
    
    for ll in 1:24
    	Avg[ll] = (sum(MessagesFired[ll,11:20])/10)/(sum(MessagesFired[ll,1:10])/10)
    end
    
    TimeSeries(Avg,NN,"Plots/TopAvg-BotAvg-PolarvsNonPolar","Number of Nodes = 20 - Power Top/Bottom vs Time")

  #  temp1 = sum(ResponsesTriggered[:,:,1],2)
  #  temp2 = sum(ResponsesTriggered[:,:,2],2)
  #  data[:,1] = temp1./temp1[1]
  #  data[:,2] = temp2./temp2[1]

  #  heatMap(ResponsesTriggered[:,:,1], "Plots/ResponsesTriggered_1-8000",1);
  #  heatMap(ResponsesTriggered[:,:,2], "Plots/ResponsesTriggered_8001-16000",2);
  #  multiHistogram(data[:,:,1], "Plots/ResponsesTriggered_GroupedBar_Lowest", "Responses Relative to Lowest - Polarized to Non-Polarized - Nodes = 20")
    
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
        BBtmp = abs.(BBtmp)
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
        BBtmp = abs.(BBtmp)
        #BB = BBtmp / norm(BBtmp)
        BB = (BB*BBtmp) / norm(BB)
    end
   return BB
end



