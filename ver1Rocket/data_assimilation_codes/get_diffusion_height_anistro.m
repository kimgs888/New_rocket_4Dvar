function [HGT_anistrx,dst_frm_src_sqd] = get_diffusion_height_anistro (temp_folder,savetime)
% This function extracts the already computed and store  distance from the
% source r^2 , the source is the emission point 
% 
    if ispc 
         slashx = '\';
    elseif isunix 
        slashx = '/';
    end
tm = num2str(savetime);
tm(tm=='.')='_';
filex = [temp_folder, slashx, 'HGX',tm,'.mat'];
HGT_anistro = load(filex);
HGT_anistrx = HGT_anistro.newdiffusion_height;
filed = [temp_folder, slashx, 'dst_frm_src_sqd',tm,'.mat'];
dst_frm_src2 = load(filed);
dst_frm_src_sqd = dst_frm_src2.dist_from_source;

clear HGT_anistro dst_frm_src2

 