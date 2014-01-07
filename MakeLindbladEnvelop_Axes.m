function MakeLindbladEnvelop_Axes (axes_Handle, a)
   
    %hold all
    
    DrawLindblad(axes_Handle, a);
    
    %set(get(axes_Handle,'Title'), 'Interpreter', 'latex');
    %title('${I} = \sqrt{{r^4}{r_0^2+r^2}^{-3/2}}, {H} = {r^2/2}{r_0^2+r^2}^{-3/2}-{r_0^2+r^2}^{-1/2}$');
    h = get(axes_Handle,'Title');
    set(h, 'Interpreter', 'tex');
    set(h, 'String', '${I} = \sqrt{{r^4}{r_0^2+r^2}^{-3/2}}, {H} = {r^2/2}{r_0^2+r^2}^{-3/2}-{r_0^2+r^2}^{-1/2}$');

    %hold off

end