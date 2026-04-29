clc; clear; close all;

%% compare between bit file and data file

filename = 'processed_bits/processed_bits12_Feb_2026_09_30_39_442.txt';
fileID = fopen(filename,'r');
formatSpec = '%s';
deECCbits = fscanf(fileID,formatSpec);
fclose(fileID);
deECCbits = dec2bin(hex2dec(reshape(deECCbits, numel(deECCbits), [])))';
deECCbits = (deECCbits(:)).';

%% parsing
payload_length = bin2dec(deECCbits(1:8));

unknown1 = bin2dec(deECCbits(9:16));

version = bin2dec(deECCbits(17:24));

sequence_number = swapbytes(uint16(bin2dec(deECCbits(25:40))));

states_info = deECCbits(41:56);

serial = deECCbits(57:184);
serial = reshape(serial, 8, [])';
serial = bin2dec(serial);
serial = char(serial); 

divider_for_degrees = 174533; % A given number that we need to divide by to get coordinate in degrees

long = swapbytes(int32(bin2dec(['0b',deECCbits(185:216),'s32'])));
long = vpa(double(long)/divider_for_degrees);

lat = swapbytes(int32(bin2dec(['0b',deECCbits(217:248),'s32'])));
lat = vpa(double(lat)/divider_for_degrees);

altitude = swapbytes(int16(bin2dec(['0b',deECCbits(249:264),'s16'])));

height =  swapbytes(int16(bin2dec(['0b',deECCbits(265:280),'s16'])));

v_north =  swapbytes(int16(bin2dec(['0b',deECCbits(281:296),'s16'])));

v_east =  swapbytes(int16(bin2dec(['0b',deECCbits(297:312),'s16'])));

v_up = swapbytes(int16(bin2dec(['0b',deECCbits(313:328),'s16'])));

unknown2 = swapbytes(int16(bin2dec(['0b',deECCbits(329:344),'s16'])));

time = datetime(1970,1,1,0,0,0,0,'InputFormat','dd-MM-uuuu'' ''HH:mm:ss');
addMS = ceil(milliseconds(swapbytes(uint64(bin2dec(deECCbits(345:408))))));
time = time + addMS;

app_lat = swapbytes(int32(bin2dec(['0b',deECCbits(409:440),'s32'])));
app_lat = vpa(double(app_lat)/divider_for_degrees);

app_long = swapbytes(int32(bin2dec(['0b',deECCbits(441:472),'s32'])));
app_long = vpa(double(app_long)/divider_for_degrees);

home_long = swapbytes(int32(bin2dec(['0b',deECCbits(473:504),'s32'])));
home_long = vpa(double(home_long)/divider_for_degrees);

home_lat = swapbytes(int32(bin2dec(['0b',deECCbits(505:536),'s32'])));
home_lat = vpa(double(home_lat)/divider_for_degrees);

device_type = bin2dec(deECCbits(537:544));

uuid_len = bin2dec(deECCbits(545:552));

uuid = deECCbits(553:712);
uuid = reshape(uuid, 8, [])';
uuid = bin2dec(uuid);
uuid = char(uuid); 

crc = uint16(bin2dec(deECCbits(713:end)));
crc = dec2hex(swapbytes(crc));

fileID = fopen('test.txt','w');
fprintf(fileID,'payload_length:%d,\n',payload_length);
fprintf(fileID,'unknown1:%d,\n',unknown1);
fprintf(fileID,'version:%d,\n',version);
fprintf(fileID,'sequence_number:%d,\n',sequence_number);
fprintf(fileID,'states_info:[%c,%c,%c,%c,%c,%c,%c,%c,%c,%c,%c,%c,%c,%c,%c,%c],\n',states_info);
fprintf(fileID,'serial:"%s",\n',serial);
fprintf(fileID,'long:%.30f,\n',long);
fprintf(fileID,'lat:%.30f,\n',lat);
fprintf(fileID,'altitude:%d,\n',altitude);
fprintf(fileID,'height:%d,\n',height);
fprintf(fileID,'v_north:%d,\n',v_north);
fprintf(fileID,'v_east:%d,\n',v_east);
fprintf(fileID,'v_up:%d,\n',v_up);
fprintf(fileID,'unknown2:%d,\n',unknown2);
fprintf(fileID,'time:%s,\n',time);
fprintf(fileID,'app_lat:%.30f,\n',app_lat);
fprintf(fileID,'app_long:%.30f,\n',app_long);
fprintf(fileID,'home_long:%.30f,\n',home_long);
fprintf(fileID,'home_lat:%.30f,\n',home_lat);
fprintf(fileID,'device_type:%d,\n',device_type);
fprintf(fileID,'uuid_len: %d,\n',uuid_len);
fprintf(fileID,'uuid: "%s",\n',uuid);
fprintf(fileID,'crc:"%s",\n',crc);
fprintf(fileID,'crc_valid:true');

fclose(fileID);