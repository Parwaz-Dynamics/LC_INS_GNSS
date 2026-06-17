import serial
import csv
import os

# ========================= CONFIGURATION =========================
PORT = 'COM7'                  # ← Change if needed
BAUDRATE = 115200              # ← Most likely value — change to 9600/38400/etc if no data appears
# ================================================================

IMU_FILE = 'IMU.csv'
GPS_FILE = 'GPS.csv'

# Headers (without the "Type" column)
IMU_HEADER = ['gps_obssec(s)', 'ax', 'ay', 'az', 'gx', 'gy', 'gz']
GPS_HEADER = ['gps_obssec(s)', 'lat', 'lng', 'alt', 'vn', 've', 'hdop', 'no. of sats']

def init_csv(filename, header):
    """Delete file if exists, then create with header"""
    if os.path.exists(filename):
        os.remove(filename)
        print(f"Deleted existing file: {filename}")
    with open(filename, 'w', newline='') as f:
        csv.writer(f).writerow(header)
    print(f"Created {filename} with header")

init_csv(IMU_FILE, IMU_HEADER)
init_csv(GPS_FILE, GPS_HEADER)

# Open serial port
try:
    ser = serial.Serial(
        port=PORT,
        baudrate=BAUDRATE,
        timeout=1,
        bytesize=serial.EIGHTBITS,
        parity=serial.PARITY_NONE,
        stopbits=serial.STOPBITS_ONE
    )
    print(f"✅ Listening on {PORT} @ {BAUDRATE} baud...")
    print("Press Ctrl+C to stop\n")

except serial.SerialException as e:
    print(f"❌ Could not open {PORT}: {e}")
    print("Check port name, device connection, and baud rate.")
    exit(1)


try:
    while True:
        line = ser.readline()

        if not line:
            continue

        try:
            # Try to decode — ignore invalid bytes
            line_str = line.decode('ascii', errors='ignore').strip()
        except Exception:
            continue

        if not line_str:
            continue

        # Split and clean fields
        fields = [f.strip() for f in line_str.split(',')]

        if len(fields) == 8 and fields[0] == 'IMU':
            data = fields[1:]  # skip 'IMU'
            with open(IMU_FILE, 'a', newline='') as f:
                csv.writer(f).writerow(data)
            # print(f"→ IMU saved   {fields[1]:>10} s")   # fixed typo

        elif len(fields) == 9 and fields[0] == 'GPS':
            data = fields[1:]  # skip 'GPS'
            with open(GPS_FILE, 'a', newline='') as f:
                csv.writer(f).writerow(data)
            # print(f"→ GPS saved   {fields[1]:>10} s")   # fixed typo

        else:
            if fields and fields[0] in ('IMU', 'GPS'):
                print(f"⚠️ Field count mismatch: {len(fields)} fields → {line_str}")
            else:
                print(f"⚠️ Skipped: {line_str}")

except KeyboardInterrupt:
    print("\n\n⛔ Stopped by user")

except Exception as e:
    print(f"\nUnexpected error: {e}")

finally:
    if 'ser' in locals() and ser.is_open:
        ser.close()
    print("Serial port closed.")