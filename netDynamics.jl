function dynamicCommunicators(NN=40,TT=365,pb=0.501,cb=1,cr=4)
    imp = [e^x for x in 1:NN]
    AA = zeros(NN,NN,TT)
    #init Data
    for ii in 1:NN
        if(rand(1)[1] < pb)
            tmp = deleteat!(collect(1:NN),ii)
            destinationN = tmp[rand(1:end)]            
            AA[ii,destinationN,1] = 1
        end

    end
    #produce the data
    #produce the Basal rate in the matrix
    for tt in 1:(TT-1)
        #Basal Loop
        for ii in 1:NN
            if(rand(1)[1] < pb)
                tmp = deleteat!(collect(1:NN),ii)
                destinationN = tmp[rand(1:end)]            
                AA[ii,destinationN,tt+1] = 1
            end

        end
        #Response Loop
        for ii in 1:NN
            msgsTo_ii = AA[:,ii]
            totalImportance = sum(msgsTo_ii .* imp)
            r_nNumerator = totalImportance
            r_nDenominator = 1 + (imp[end] * sum(msgsTo_ii))
            r_next = r_nNumerator / r_nDenominator
            if(r_next > rand(1)[1])#if so generate cr links
                println("response")
                tmp = deleteat!(collect(1:NN),ii)
                destinationNodes = tmp[randperm(length(tmp))[1:cr]]
                AA[ii,destinationNodes,tt+1] = 1

            end
        end
                
    end
    print(AA)

end
