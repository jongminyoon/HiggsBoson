# const v = 0.326714
# const mw = 0.103773
# const mz = 0.119549
# const mt = 0.191556
# const mh = 0.13573 

# const v = 0.270619
# const mw = 0.085902 
# const mz = 0.098971  
# const mt = 0.160562
# const mb = 0.00292719
# const mh = 0.135871 

# const g = 2*mw/v; 
# const gy = sqrt((2*mz/v)^2-g^2);


# const gy = 0.364214    
# const g = 0.633234    
# const yt = 0.814796    
# # const vev  = 272.15/1000
# const vev  = 246/1000
# const mw = g/2*vev
# const mz = sqrt(gy^2 + g^2)/2*vev
# const mt = yt/sqrt(2)*vev
# # const mh = 136.64/1000
# const mh = 125.09/1000





const vev  = 246/1000

const yt = 0.843514
const yb = 0.0152661
const ytau = 0.0106599

const mt = yt/sqrt(2)*vev
const mb = yb/sqrt(2)*vev
const mtau = ytau/sqrt(2)*vev

const mw = 80.384/1000
const mz = 91.1876/1000
const g = 2*mw/vev
const gy = sqrt((2*mz/vev)^2 - g^2)    

const mh = 125.09/1000




const zratio = 1/100; 
const kpir = -log(zratio);