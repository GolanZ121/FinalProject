clc; clear; close all;

%% tests
deECCbits = ones(1,1000);
payload_length = bin2dec(int2str(deECCbits(1:8)));
unknown1 = bin2dec(int2str(deECCbits(9:16)));
version = bin2dec(int2str(deECCbits(17:24)));
sequence_number = bin2dec(int2str(deECCbits(25:40)));
states_info = deECCbits(41:56);
serial = deECCbits(57:184);
serial = reshape(serial, 8, [])';
serial = bin2dec(int2str(serial));
serial = char(serial); 
long = bin2dec(int2str(deECCbits(185:216)));
lat = bin2dec(int2str(deECCbits(217:248)));
altitude = bin2dec(int2str(deECCbits(249:264)));
height = bin2dec(int2str(deECCbits(265:280)));
v_north = bin2dec(int2str(deECCbits(281:296)));
v_east = bin2dec(int2str(deECCbits(297:312)));
v_up = bin2dec(int2str(deECCbits(313:328)));
unknown2 = bin2dec(int2str(deECCbits(329:344)));
% time = %COMPLETE LATER!!!
app_lat = bin2dec(int2str(deECCbits(409:440)));
app_long = bin2dec(int2str(deECCbits(441:472)));
home_long = bin2dec(int2str(deECCbits(473:504)));
home_lat = bin2dec(int2str(deECCbits(505:536)));
device_type = bin2dec(int2str(deECCbits(537:544)));
uuid_len = bin2dec(int2str(deECCbits(545:552)));
uuid = deECCbits(57:184);
uuid = reshape(uuid, 8, [])';
uuid = bin2dec(int2str(uuid));
uuid = char(uuid); 

fileID = fopen('test.txt','w');
fprintf(fileID,'payload_length:%d,\n',payload_length);
fprintf(fileID,'unknown1:%d,\n',unknown1);
fprintf(fileID,'version:%d,\n',version);
fprintf(fileID,'sequence_number:%d,\n',sequence_number);
fprintf(fileID,'states_info:[%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d],\n',states_info);
fprintf(fileID,'serial:"%s",\n',serial);
fprintf(fileID,'long:%d,\n',long);
fprintf(fileID,'lat:%d,\n',lat);
fprintf(fileID,'altitude:%d,\n',altitude);
fprintf(fileID,'height:%d,\n',height);
fprintf(fileID,'v_north:%d,\n',v_north);
fprintf(fileID,'v_east:%d,\n',v_east);
fprintf(fileID,'v_up:%d,\n',v_up);
fprintf(fileID,'unknown2:%d,\n',unknown2);
fprintf(fileID,'time:%d-%d-%d,\n',unknown2);%COMPLETE LATER!!!
fprintf(fileID,'app_lat:%d,\n',app_lat);
fprintf(fileID,'app_long:%d,\n',app_long);
fprintf(fileID,'home_long:%d,\n',home_long);
fprintf(fileID,'home_lat:%d,\n',home_lat);
fprintf(fileID,'device_type:%d,\n',device_type);
fprintf(fileID,'uuid_len:%d,\n',uuid_len);
fprintf(fileID,'uuid:"%s",\n',serial);

fclose(fileID);
 
% packet = [
%     "time":"12-Feb-2026 07:31:30","crc":"38F4","crc_valid":true]