function dynamicCommunicators(NN,tt)
    pb = 0.25
    cb = 1
    cr = 4
    AA = zeros(NN,NN,tt)
    #produce the data
    #produce the Basal rate in the matrix
    for ii in 1:NN
        if(rand(1)[1] < pb)#produce Basal
            tmp = deleteat!(collect(1:NN),ii)
            destinationN = tmp[rand(1:end)]            
            AA[ii,destinationN,1] = 1
        end
        
    end
    
    print(AA)

end
