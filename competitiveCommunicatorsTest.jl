using GR
#pyplot()
function main(NN=5,TT=3,pb=0.6,cb=1,cr=4)
    
    AA = dynamicCommunicators(NN,TT,pb,cb,cr)
 #  println(AA, " dynamicCommunicators worked")
    degOut = totalOut(AA,NN)
 #  println(degOut, " totalOut worked")
    Cbroad = dynamicCentrality(AA,NN,TT)
    println(Cbroad')
 #  scatterPlot(AA,size(AA)[1],degOut,Cbroad)

   # AA2 = ringNet()
   # degOut = totalOut(AA2,21)
   # Cbroad = dynamicCentrality(AA2,21,3)
   # scatterPlot(AA2,size(AA2)[1],degOut,Cbroad)
end


function scatterPlot(AA,NN,degOut,Cbroad)
    labels = [string(x) for x in 1:NN]
    GR.scatter(degOut,Cbroad,markersize=1,series_annotations =labels,grid=false)
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


function dynamicCommunicators(NN=40,TT=365,pb=0.01,cb=1,cr=4)
    imp = ([e^x for x in 1:NN] .* ones(NN,NN))';
    println("imp = ", imp)

    AA = zeros(NN,NN,TT)
    #init Data
    for ii in 1:NN
        if(rand(1)[1] < pb)
            tmp = deleteat!(collect(1:NN),ii)
            destinationN = tmp[rand(1:end)]            
            AA[ii,destinationN,1] = 1
        end

    end

    # Change importance, produce basal, and generate response
    for tt in 1:(TT-1)
	#alter importance based on current matrix state before generating new
 	for ii in 1:(NN-1)
	    for jj in 1:(NN-1)
		if(AA[ii,jj,TT] == 0) #&& AA[jj,ii,TT] == 0)
			imp[ii,jj]=0;
		#	imp[jj,ii]*=.5;
		end
	    end
	end
        #Basal Loop
        for ii in 1:NN
            if(rand(1)[1] < pb)
                #print("basal")
                tmp = deleteat!(collect(1:NN),ii)
                destinationN = tmp[rand(1:end)]            
                AA[ii,destinationN,tt+1] = 1
            end

        end
        #Response Loop
        for ii in 1:NN
            msgsTo_ii = AA[:,ii,tt]
            
            totalImportance = sum(msgsTo_ii .* imp)
            
            r_nNumerator = totalImportance
            r_nDenominator = 1 + (e^NN * sum(msgsTo_ii))
            r_next = r_nNumerator / r_nDenominator
            #println(r_next)
            if(r_next >  rand(1)[1])#if so generate cr links
                #print("response=");print(find(totalImportance));println(":")
                # println(msgsTo_ii)
                tmp = deleteat!(collect(1:NN),ii)
                destinationNodes = tmp[randperm(length(tmp))[1:cr]]
                AA[ii,destinationNodes,tt+1] = 1

            end
        end
        println(AA)       
    end
    println(imp)
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
