using Plots
function dynamicCommunicators(NN=40,TT=365,pb=0.1,cb=1,cr=4)
    imp = [x^4 for x in 1:NN]
    println(imp)
    AA = zeros(NN,NN,TT)
    #init Data
    for ii in 1:NN
        if(rand(1)[1] < pb)
            tmp = deleteat!(collect(1:NN),ii)
            destinationN = tmp[rand(1:end)]            
            AA[ii,destinationN,1] = 1
        end

    end
    #produce the Basal rate in the matrix
    for tt in 1:(TT-1)
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
            r_nDenominator = 1 + (findmax(imp)[1] * sum(msgsTo_ii))
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
                
    end
    #print(AA)
    #now get the total out degree of the nodes
    degOut = totalOut(AA,NN)
    println("degOut")
    print(degOut)
    Cbroad = dynamicCentrality(AA,NN,TT)
    labels = [string(x) for x in 1:NN]
    
    scatter(degOut,Cbroad,markersize=1,series_annotations =labels,grid=false)
    
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
    for kk = 1:TT#:-1:1
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
    #println(Cbroad)
    #println(0.9*(1/findmax(daynorm)[1]))
    return Cbroad
end


#for k = size(Atensor,3):-1:1,
#    Aday = Atensor(:,:,k);
#    B = (eye(size(Atensor,1),size(Atensor,1)) - alpha*Aday)\B;
#    B = abs(B);
#    B = B/norm(B);
#end
