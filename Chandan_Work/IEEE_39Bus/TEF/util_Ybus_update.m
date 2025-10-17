% util_Ybus_update.m
% Fast update of Ybus after a simple topology change (line removal or bus fault)
function Ynew = util_Ybus_update(Ybus, faultInfo)
Ynew = Ybus;
% Example placeholder: if faultInfo.type indicates line outage, zero corresponding elements
% TODO: parse faultInfo.location and update Y accordingly
end
