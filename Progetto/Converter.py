import meshio

input_path_Geo_0 = "/Users/gingerale/Documents/C++/Docker/Data/ProgettoPCS2025/Progetto/Output/Geodetico/GeodeticoCell0Ds.inp"
output_path_Geo_0 = "/Users/gingerale/Downloads/GeodeticoCell0Ds.vtu"
input_path_Geo_1 = "/Users/gingerale/Documents/C++/Docker/Data/ProgettoPCS2025/Progetto/Output/Geodetico/GeodeticoCell1Ds.inp"
output_path_Geo_1 = "/Users/gingerale/Downloads/GeodeticoCell1Ds.vtu"
input_path_Ori_0 = "/Users/gingerale/Documents/C++/Docker/Data/ProgettoPCS2025/Progetto/Output/Geodetico/OriginaleCell0Ds.inp"
output_path_Ori_0 = "/Users/gingerale/Downloads/OriginaleCell0Ds.vtu"
input_path_Ori_1 = "/Users/gingerale/Documents/C++/Docker/Data/ProgettoPCS2025/Progetto/Output/Geodetico/OriginaleCell1Ds.inp"
output_path_Ori_1 = "/Users/gingerale/Downloads/OriginaleCell1Ds.vtu"


print("Lettura della mesh in corso...")
mesh = meshio.read(input_path_Geo_0, file_format="avsucd")
print("Salvataggio nel formato VTU...")
meshio.write(output_path_Geo_0, mesh)

print("Lettura della mesh in corso...")
mesh = meshio.read(input_path_Geo_1, file_format="avsucd")
print("Salvataggio nel formato VTU...")
meshio.write(output_path_Geo_1, mesh)

print("Lettura della mesh in corso...")
mesh = meshio.read(input_path_Ori_0, file_format="avsucd")
print("Salvataggio nel formato VTU...")
meshio.write(output_path_Ori_0, mesh)

print("Lettura della mesh in corso...")
mesh = meshio.read(input_path_Ori_1, file_format="avsucd")
print("Salvataggio nel formato VTU...")
meshio.write(output_path_Ori_1, mesh)

risposta = input("Goldberg [y/n]: ").strip().lower()

if risposta == 'y':

    input_path_Gol_0 = "/Users/gingerale/Documents/C++/Docker/Data/ProgettoPCS2025/Progetto/Output/Goldberg/GoldbergCell0Ds.inp"
    output_path_Gol_0 = "/Users/gingerale/Downloads/GoldbergCell0Ds.vtu"
    input_path_Gol_1 = "/Users/gingerale/Documents/C++/Docker/Data/ProgettoPCS2025/Progetto/Output/Goldberg/GoldbergCell1Ds.inp"
    output_path_Gol_1 = "/Users/gingerale/Downloads/GoldbergCell1Ds.vtu"
    input_path_Dua_0 = "/Users/gingerale/Documents/C++/Docker/Data/ProgettoPCS2025/Progetto/Output/Goldberg/DualeCell0Ds.inp"
    output_path_Dua_0 = "/Users/gingerale/Downloads/DualeCell0Ds.vtu"
    input_path_Dua_1 = "/Users/gingerale/Documents/C++/Docker/Data/ProgettoPCS2025/Progetto/Output/Goldberg/DualeCell1Ds.inp"
    output_path_Dua_1 = "/Users/gingerale/Downloads/DualeCell1Ds.vtu"

    print("Lettura della mesh in corso...")
    mesh = meshio.read(input_path_Gol_0, file_format="avsucd")
    print("Salvataggio nel formato VTU...")
    meshio.write(output_path_Gol_0, mesh)

    print("Lettura della mesh in corso...")
    mesh = meshio.read(input_path_Gol_1, file_format="avsucd")
    print("Salvataggio nel formato VTU...")
    meshio.write(output_path_Gol_1, mesh)

    print("Lettura della mesh in corso...")
    mesh = meshio.read(input_path_Dua_0, file_format="avsucd")
    print("Salvataggio nel formato VTU...")
    meshio.write(output_path_Dua_0, mesh)

    print("Lettura della mesh in corso...")
    mesh = meshio.read(input_path_Dua_1, file_format="avsucd")
    print("Salvataggio nel formato VTU...")
    meshio.write(output_path_Dua_1, mesh)

elif risposta == 'n':
    print("")

else:
    print("Errore: Input non riconosciuto. Per favore esegui di nuovo lo script e inserisci solo 'y' o 'n'.")


print("Conversione completata!")