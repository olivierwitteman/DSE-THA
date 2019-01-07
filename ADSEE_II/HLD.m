function [Delta_Clmax_HLD] = HLD(LE, TE , c_ac_c)
Delta_Clmax_TE=0;
Delta_Clmax_LE=0.;
disp('HLD')
if strcmp(TE,'plain_split')
    Delta_Clmax_TE=0.9;
end
if strcmp(TE,'Slotted')
    Delta_Clmax_TE=1.3;
end 
if strcmp(TE,'fowler')==1.
    Delta_Clmax_TE=1.3*c_ac_c;
end
if strcmp(TE,'double_slotted')==1.
    Delta_Clmax_TE=1.6*c_ac_c;
end
if strcmp(TE, 'triple_slotted')==1.
    Delta_Clmax_TE=1.9*c_ac_c;
end
if strcmp(LE,'fixed_slot')==1.
    Delta_Clmax_LE=0.2;
end
if strcmp(LE,'LE_flap')==1.
    Delta_Clmax_LE=0.3;
end
if strcmp(LE,'kruger_flap')==1.
    Delta_Clmax_LE=0.3;
end
if strcmp(LE,'Slat')==1.
    Delta_Clmax_LE=0.4*c_ac_c;
end

Delta_Clmax_HLD= Delta_Clmax_LE +Delta_Clmax_TE
end