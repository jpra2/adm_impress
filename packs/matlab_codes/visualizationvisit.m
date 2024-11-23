function visualizationvisit(S_old,vpi_old,p,time2)

if vpi_old >= time2
         step = 100*time2;
        creat_vtk_file(p,S_old,step);
        %creat_vtk_file_quad(coord,elem,p,S_old,step)
         time2 = time2 + 0.01;
end


end