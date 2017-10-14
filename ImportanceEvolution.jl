using Gadfly

function main(NN=40,TT=365,pb=0.1,cb=1,cr=4)
    
    AA = dynamicCommunicators(NN,TT,pb,cb,cr)
    degOut = totalOut(AA,NN)
    Cbroad = dynamicCentrality(AA,NN,TT)
    scatterPlot(AA,size(AA)[1],degOut,Cbroad)

   # AA2 = ringNet()
   # degOut = totalOut(AA2,21)
   # Cbroad = dynamicCentrality(AA2,21,3)
   # scatterPlot(AA2,size(AA2)[1],degOut,Cbroad)
end


function scatterPlot(AA,NN,degOut,Cbroad)
    labels = [string(x) for x in 1:NN]
    plot(x = degOut,y = Cbroad,label = labels, Geom.label)
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


function dynamicCommunicators(NN=40,TT=50,pb=0.01,cb=1,cr=4)
    imp = [x^4 for x in 1:NN]
 #  imp[NN] = 100;
    initial_imp = imp 
    println(imp)
    AA = zeros(NN,NN,TT)
    #INIT DATA>>
    for ii in 1:NN
        if(rand(1)[1] <= pb)
            tmp = deleteat!(collect(1:NN),ii)
            destinationN = tmp[rand(1:end)]            
            AA[ii,destinationN,1] = 1
        end
	println(AA[:,:,1])
    end
    #DATA LOOP: 1st put basal edges into <present> / 2nd put response from past to now
    for tt in 2:(TT)
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
                totalImportance = sum( msgsTo_ii .* imp  )
                r_nNumerator = totalImportance
                r_nDenominator = 1 + (findmax(imp)[1] * sum(AA[:,ii,tt-1]))
                println("prob of response = ", r_nNumerator / r_nDenominator )
                r_next = r_nNumerator / r_nDenominator

                if(r_next >  rand(1)[1])#if so generate cr links

                    tmp = deleteat!(collect(1:NN),ii)
                    destinationNodes = tmp[randperm(length(tmp))[1:cr]]
                    AA[ii,destinationNodes,tt] = 1#FIRE NOW
		    imp = (imp += .1*msgsTo_ii);  #Increase successful nodes importance by .1

               end
            end
        end
       #  println(AA[:,:,tt], "\n")
	# println("Imp = ", imp, "\n")       
    end
    for ii in 1:NN
	println(ii, " initial_imp = ",initial_imp[ii], "   final_imp = ", imp[ii], " Responses generated = ", 10*(imp[ii]-initial_imp[ii]), "\n\n")
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

