#projection sketches

for nx in 1:Nx
    for ny in 1:Nx
        for nz in 1:Nx
#            #kx = sin(2*pi*nvalues[nx]/Nx) # for symmetric derivatives
#            #ky = sin(2*pi*nvalues[ny]/Nx) # for symmetric derivatives
#            #kz = sin(2*pi*nvalues[nz]/Nx) # for symmetric derivatives
#            #kx = 2*sin(pi*nvalues[nx]/Nx) # another derivative
#            #ky = 2*sin(pi*nvalues[ny]/Nx) # another derivative
#            #kz = 2*sin(pi*nvalues[nz]/Nx) # another derivative
#            #kx = -im*( exp(im*2*pi*nvalues[nx]/Nx)-1 )
#            #ky = -im*( exp(im*2*pi*nvalues[ny]/Nx)-1 )
#            #kz = -im*( exp(im*2*pi*nvalues[nz]/Nx)-1 )
#            #k2 = real( conj(kx)*kx + conj(ky)*ky + conj(kz)*kz )
            k = -im*[ exp(-im*2*pi*nvalues[nx]/Nx)-1 , exp(-im*2*pi*nvalues[ny]/Nx)-1 , exp(-im*2*pi*nvalues[nz]/Nx)-1]
            k2 = real( sum( abs.(k).^2 ))
            DLk[nx,ny,nz] = 0
            DTk[nx,ny,nz] = 0
#
           if k2 > 0
               for i in 1:3
                   for j in 1:3
                       DLk[nx,ny,nz] += k[i]*Ai1k[j][nx,ny,nz] * conj( k[j]*Ai1k[i][nx,ny,nz] ) / k2
                   end
                   DTk[nx,ny,nz] += Ai1k[i][nx,ny,nz] * conj(Ai1k[i][nx,ny,nz])
               end
               DTk[nx,ny,nz] -= DLk[nx,ny,nz]
           end
#
#            if k2 > 0 
#                #DLk[nx,ny,nz] = ( (kx)*conj(kx)*Dk[1,1][nx,ny,nz] + (ky)*conj(ky)*Dk[2,2][nx,ny,nz] + (kz)*conj(kz)*Dk[3,3][nx,ny,nz]
#                #                + (kx)*conj(ky)*Dk[1,2][nx,ny,nz] + (ky)*conj(kx)*Dk[2,1][nx,ny,nz] 
#                #                + (kx)*conj(kz)*Dk[1,3][nx,ny,nz] + (kz)*conj(kx)*Dk[3,1][nx,ny,nz] 
#                #                + (ky)*conj(kz)*Dk[2,3][nx,ny,nz] + (kz)*conj(ky)*Dk[3,2][nx,ny,nz] ) / k2
#                #DLk[nx,ny,nz] = ( k[1]*Ai1k[1]*conj(k[1]*Ai1k[1]) + k[2]*Ai1k[]*conj(k[2]*Ai1k[]) + k[3]*Aik1[]*conj(kz*Aik1[])
#                #                + k[1]*Ai1k[1]*conj(k[2]*Ai1k[2]) + k[2]*Ai1k[]*conj(k[1]*Ai1k[]) 
#                #                + k[1]*Ai1k[1]*conj(k[3]*Ai1k[3]) + k[3]*Ai1k[]*conj(k[1]*Ai1k[]) 
#                #                + k[2]*Ai1k[2]*conj(k[3]*Ai1k[3]) + k[3]*Ai1k[]*conj(k[2]*Ai1k[]) ) / k2
#                #DLk[nx,ny,nz] = ( (kx)*Dk[1,1][nx,ny,nz] + (ky)*Dk[2,2][nx,ny,nz] + (kz)*Dk[3,3][nx,ny,nz]
#                #                + (kx)*Dk[1,2][nx,ny,nz] + (ky)*Dk[2,1][nx,ny,nz] 
#                #                + (kx)*Dk[1,3][nx,ny,nz] + (kz)*Dk[3,1][nx,ny,nz] 
#                #                + (ky)*Dk[2,3][nx,ny,nz] + (kz)*Dk[3,2][nx,ny,nz] ) / k2
#                #
#                #DTk[nx,ny,nz] = ( Dk[1,1][nx,ny,nz] + Dk[2,2][nx,ny,nz] + Dk[3,3][nx,ny,nz] )
#                #DTk[nx,ny,nz] -= (conj(kx)*kx*Dk[1,1][nx,ny,nz] + conj(ky)*ky*Dk[2,2][nx,ny,nz] + conj(kz)*kz*Dk[3,3][nx,ny,nz]
#                #                + conj(kx)*ky*(Dk[1,2][nx,ny,nz]+Dk[2,1][nx,ny,nz]) 
#                #                + conj(kx)*kz*(Dk[1,3][nx,ny,nz]+Dk[3,1][nx,ny,nz]) 
#                #                + conj(ky)*kz*(Dk[2,3][nx,ny,nz]+Dk[3,2][nx,ny,nz]) ) / k2
#            end
       end
   end
end
#@show DLk[123]
DLk[1] = ( Dk[1,1][1] + Dk[2,2][1] + Dk[3,3][1] + Dk[1,2][1] + Dk[2,1][1] + Dk[1,3][1] + Dk[3,1][1]  + Dk[2,3][1] + Dk[3,2][1] )
DTk[1] = ( Dk[1,1][1] + Dk[2,2][1] + Dk[3,3][1] - DLk[1] ) 
#DTk *= 1/(disc.sdim-1) # divide by transverse projection normalisation P_TP_T = sdim - 1

#@show sum(abs.(DLk))
#@show sum(abs.(DTk))

# for nx in 1:Nx
#     for ny in 1:Nx
#         for nz in 1:Nx
#             kx = 2*exp(-im*pi*nx)*sin(pi*nx/Nx)
#             ky = 2*exp(-im*pi*ny)*sin(pi*ny/Nx)
#             kz = 2*exp(-im*pi*nz)*sin(pi*nz/Nx)
#             k2 = real( conj(kx)*kx + conj(ky)*ky + conj(kz)*kz )
            
#             if k2 > 0
#                 DLk[nx,ny,nz] = ( conj(kx)*kx*Dk[1,1][nx,ny,nz] + conj(ky)*ky*Dk[2,2][nx,ny,nz] + conj(kz)*kz*Dk[3,3][nx,ny,nz]
#                                 + conj(kx)*ky*(Dk[1,2][nx,ny,nz]+Dk[2,1][nx,ny,nz]) 
#                                 + conj(kx)*kz*(Dk[1,3][nx,ny,nz]+Dk[3,1][nx,ny,nz]) 
#                                 + conj(ky)*kz*(Dk[2,3][nx,ny,nz]+Dk[3,2][nx,ny,nz]) ) / k2
#                 DTmatk[1,1][nx,ny,nz] = Dk[1,1][nx,ny,nz] - conj(kx)*kx * DLk[nx,ny,nz] / k2
#                 DTmatk[1,2][nx,ny,nz] = Dk[1,2][nx,ny,nz] - conj(kx)*ky * DLk[nx,ny,nz] / k2
#                 DTmatk[2,1][nx,ny,nz] = Dk[2,1][nx,ny,nz] - conj(ky)*kx * DLk[nx,ny,nz] / k2
#                 DTmatk[2,2][nx,ny,nz] = Dk[2,2][nx,ny,nz] - conj(ky)*ky * DLk[nx,ny,nz] / k2
#                 #
#                 DTmatk[1,3][nx,ny,nz] = Dk[1,3][nx,ny,nz] - conj(kx)*kz * DLk[nx,ny,nz] / k2
#                 DTmatk[3,1][nx,ny,nz] = Dk[3,1][nx,ny,nz] - conj(kz)*kx * DLk[nx,ny,nz] / k2
#                 DTmatk[2,3][nx,ny,nz] = Dk[2,3][nx,ny,nz] - conj(ky)*kz * DLk[nx,ny,nz] / k2
#                 DTmatk[3,2][nx,ny,nz] = Dk[3,2][nx,ny,nz] - conj(kz)*ky * DLk[nx,ny,nz] / k2
#                 DTmatk[3,3][nx,ny,nz] = Dk[3,3][nx,ny,nz] - conj(kz)*kz * DLk[nx,ny,nz] / k2
#                 #E2long[nx,ny,nz] = (  conj(kx)*kx*E2k[1,1][nx,ny,nz] + conj(ky)*ky*E2k[2,2][nx,ny,nz] + conj(kz)*kz*E2k[3,3][nx,ny,nz]
#                 #                    + conj(kx)*ky*(E2k[1,2][nx,ny,nz]+E2k[2,1][nx,ny,nz]) 
#                 #                    + conj(kx)*kz*(E2k[1,3][nx,ny,nz]+E2k[3,1][nx,ny,nz]) 
#                 #                    + conj(ky)*kz*(E2k[2,3][nx,ny,nz]+E2k[3,2][nx,ny,nz]) ) / k2
#                 #E2trans[1,1][nx,ny,nz] = E2k[1,1][nx,ny,nz] - conj(kx)*kx * E2long[nx,ny,nz] / k2
#                 #E2trans[1,2][nx,ny,nz] = E2k[1,2][nx,ny,nz] - conj(kx)*ky * E2long[nx,ny,nz] / k2
#                 #E2trans[2,1][nx,ny,nz] = E2k[2,1][nx,ny,nz] - conj(ky)*kx * E2long[nx,ny,nz] / k2
#                 #E2trans[2,2][nx,ny,nz] = E2k[2,2][nx,ny,nz] - conj(ky)*ky * E2long[nx,ny,nz] / k2
#             end
             
#         end
#     end
# end

# for nx in 1:Nx
#     for ny in 1:Nx
#         for nz in 1:Nx
#             kx = 2*exp(-im*pi*nx)*sin(pi*nx/Nx)
#             ky = 2*exp(-im*pi*ny)*sin(pi*ny/Nx)
#             kz = 2*exp(-im*pi*nz)*sin(pi*nz/Nx)
#             k2 = real( conj(kx)*kx + conj(ky)*ky + conj(kz)*kz )
            
#             if k2 > 0
#                 DTk[nx,ny,nz] = DTmatk[1,1][nx,ny,nz]+DTmatk[2,2][nx,ny,nz]+DTmatk[3,3][nx,ny,nz]
#                 DTk[nx,ny,nz] -= ( conj(kx)*kx* DTmatk[1,1][nx,ny,nz] + conj(ky)*ky*DTmatk[2,2][nx,ny,nz] + conj(kz)*kz*DTmatk[3,3][nx,ny,nz]
#                                  + conj(kx)*ky*(DTmatk[1,2][nx,ny,nz]+DTmatk[2,1][nx,ny,nz]) 
#                                  + conj(kx)*kz*(DTmatk[1,3][nx,ny,nz]+DTmatk[3,1][nx,ny,nz]) 
#                                  + conj(ky)*kz*(DTmatk[2,3][nx,ny,nz]+DTmatk[3,2][nx,ny,nz]) ) / k2
#             end
#         end
#     end
# end
#@show sum(imag(DLk))
#@show sum(imag(DTk))