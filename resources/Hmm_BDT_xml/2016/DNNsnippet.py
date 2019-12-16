
### suppose that data is a numpy array produced by tree2array via root-numpy
def make_array(data, withMass, flags):
    arr = []
    arr.append((data[:]['hmass']-115.0)/20.0)
    arr.append( data[:]['nsajets2']/50.0)
    arr.append( data[:]['sajetsht2']/2000.0)
    arr.append( data[:]['hmmpt']/6500.0)
    arr.append((data[:]['hmmrap']+5.0)/10.0)
    arr.append((data[:]['hmmthetacs']+1.0)/2.0)
    arr.append((data[:]['hmmphics']+3.5)/7.0)
    arr.append( data[:]['j1pt']/6500.0)
    arr.append( data[:]['j2pt']/6500.0)
    arr.append( data[:]['j1eta']+5.0/10.0)
    arr.append( data[:]['detajj']/10.0)
    arr.append( data[:]['dphijj']/3.5)
    arr.append( data[:]['mjj']/6500.0)
    arr.append((data[:]['zepen']+15.0)/30.0)
    arr.append( data[:]['detammj']/10.0)
    arr.append( data[:]['dphimmj']/3.5)
    arr.append( data[:]['m1pt/hmass']/50.0)
    arr.append( data[:]['m2pt/hmass']/50.0)
    arr.append((data[:]['m1eta']+2.5)/5.0)
    arr.append((data[:]['m2eta']+2.5)/5.0)
    arr.append( data[:]['njets']/50.0)
    arr.append( data[:]['ptbalance'])
    arr.append( data[:]['ptcentrality']/200/)
    arr = np.asarray(arr).transpose()

    array = []
    for i in range(data.shape[0]): ### in case you want to check on testing events
        if flags[data[i]['event'] % len(flags)] > 0 : array.append(list(arr[i]))

    return np.asarray(array)
