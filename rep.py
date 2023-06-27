a={'Thejournal of Physical Chemistry Letters': 'The Journal of Physical Chemistry Letters', 'Inorganic chemistry': 'Inorganic Chemistry', 'Modelling and simulation in materials science and engineering': 'Modelling and Simulation in Materials Science and Engineering', 'Optics letters': 'Optics Letters', 'Thejournal of Physical Chemistry C': 'The Journal of Physical Chemistry C', 'RSC advances': 'RSC Advances', 'ACS Applied Materials \\& Interfaces': 'ACS Applied Materials and Interfaces', 'Physical review letters': 'Physical Review Letters', 'Journal of alloys and compounds': 'Journal of Alloys and Compounds', 'Journal of applied physics': 'Journal of Applied Physics', 'Optical materials': 'Optical Materials', 'Solid state communications': 'Solid State Communications', 'physica status solidi (b)': 'Physica Status Solidi (B)', 'Newjournal of Physics': 'New Journal of Physics', 'Chemical Engineeringjournal': 'Chemical Engineering Journal', 'Thejournal of Chemical Physics': 'The Journal of Chemical Physics', 'Nature communications': 'Nature Communications'}

text=""
flag=0
with open("ref.bib","r") as f:
    for i in f.readlines():
        print(i)
        if "journal" in i:
            for j in a:
                if j in i:
                    text+=i.replace(j,a[j])
                    flag=1
                    break
            if flag==0:
                text += i
            else:
                flag=0

        else:
            text+=i
with open("lala.bib","w") as f:
    f.write(text)



