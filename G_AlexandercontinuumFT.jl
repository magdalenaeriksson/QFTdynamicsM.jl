using FFTW

#discrete sine II

N=12
r=rand(N)
@show r

# no plan
@show FFTW.r2r(r,FFTW.RODFT10)

# with plan
#A = plan_rfft(r,FFTW.RODFT00)
B       = FFTW.plan_r2r(r, FFTW.RODFT10);   
Binv    = FFTW.plan_r2r(r, FFTW.RODFT01);  # inverse
@show B*r;
@show Binv*(B*r);

# Normalisation
@show Binv*(B*r)/(N*2) - r;

####
nrmomenta = 16
F = rand(nrmomenta) 

function calcFx(F, nrmomenta)
    isofactor = [i for i in 1:nrmomenta]
    normfactor = (1/(4*nrmomenta^2)) * [1/(i-0.5) for i in 1:nrmomenta]
    result = FFTW.r2r(isofactor .* F,FFTW.RODFT01)
    return normfactor .* result
end

function calcF(Fx, nrmomenta)
    isofactor = [(i-0.5) for i in 1:nrmomenta]
    normfactor = 2*nrmomenta * [1/i for i in 1:nrmomenta]
    result = FFTW.r2r(isofactor .* Fx,FFTW.RODFT10)
    return normfactor .* result
end

@show F
Fx = calcFx(F, nrmomenta)
@show Fx
F2 = calcF(Fx, nrmomenta)
@show F2
@show F2./F

