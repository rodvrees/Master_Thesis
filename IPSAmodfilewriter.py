with open("unimodptms.txt") as fn:
    with open("IPSA_MODfile.csv", "w") as ofl:
        ofl.write("Modification Name,Modification Mass\n")
        for line in fn.readlines():
            cutline = line[4:]
            lijst = cutline.split(",")
            AA = lijst[3]
            if AA == "N-term":
                AA = "N-TERM"
            if AA == "C-term":
                AA = "C-TERM"
            Delta_mass = lijst[1]
            mod = lijst[0]
            mod = mod.replace(":","|")
            mod_formatted = mod + "[" + AA + "]"
            ofl.write("{},{}\n".format(mod_formatted, Delta_mass))

ofl.close()
fn.close()


