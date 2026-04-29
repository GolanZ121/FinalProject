
weights = [8 4 2 1];
decimalValues = reshape(data, 4, []).' * weights';

% Convert the resulting decimals to hex characters
H = dec2hex(decimalValues).'; % Returns 'E3'

% hex_data = large_hex_data_to_bits(data);
hex_known = '5810024D00071D334E33424A393430313230344A44000000000000000000003E00000000000700FAFF2BDC371CC3509C010000CB225500CCCE5C0000000000000000003A133139383835373531313131373735323332303000F438';
bits_known = large_hex_data_to_bits(hex_known);

function raw_bits = large_hex_data_to_bits(raw_hex)
% Map each hex digit to its 4 bits binary form
hex_map = dec2bin(0:15, 4) - '0';

% Convert hex characters to indices (0-15)
idx = hex2dec(raw_hex(:)) + 1;

% Retrieve bits and reshape into a single row
raw_bits = reshape(hex_map(idx, :).', 1, []);
end
