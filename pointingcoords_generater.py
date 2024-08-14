import pandas as pd

# tele = Pointing(182.13737, -23.78930, '2028-06-01T00:00:00', '2028-06-30T23:59:59', '1h')
#
# dict = {'RA': tele.pointing_RAs, 'DEC': tele.pointing_DECs}
# df = pd.DataFrame(dict)
# df.to_csv('pointingcoords_2028.csv')

df = pd.read_csv("pointingcoords_2028.csv")
print(df.columns)