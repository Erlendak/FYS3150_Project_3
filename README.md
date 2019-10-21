# FYS3150_Project_3
The 3rd project, in FYS3150

I dette prosjektet tar vi for oss numerisk integrasjon ved hjelp av fire forskjellige metoder. Dette er Gauss Quadrature med Legendre polynomer og Laguerre polynomer og Monte Carlo med brute force og med importance sampling. Vi setter først opp funksjonen som vi skal integrere, og justerer den slik at vi kan integrere den numerisk. Dette blir et seks-dimensjonalt integral, og det er her hele poenget i prosjektet ligger. Vi vil vite hvilken metode som er best til å håndtere et slikt seks-dimensjonalt integral. Vi kan da sette opp metodene våre, vi setter først opp Gauss-Legendre ganske direkte og regner integralet. Ut fra denne metoden kan vi endrekoordiantene, polynomene og justerer funksjonen slik at vi får Gauss-Laguerre metoden, som vi også kan bruke for å løse integralet.

Vi setter så opp brute force Monte Carlo metoden direkte med uniform distribusjon i kartesiske koordinater og bruker dette til å finne ut av integralet. Og til slutt justerer vi på Monte Carlo metoden ved å bruke sfæriske koordinater, legge til en eksponential distribusjon og ved å justere funksjonen vår siden vi kan absorbere ledd inn i den eksponensiale distribusjonen. Vi parallelliserer Monte Carlo metodene for å gjøre koden mer effektiv. Vi setter også inn compiler flags for å prøve å gjøre koden enda mer effektiv.

Når vi da får ut alle resultatene som vi trenger i form av tabeller og plott kan vi sammenlikne effektiviteten og presisjonen til de forskjellige metodene med hverandre for å se hvilken metode som er mest egnet til å løse det seks-dimensjonale integralet. Vi finner ut at Importance Sampling er den mest effektive of presise metoden av de fire, og at med parallellisering så er veldig rask. Vi fant også ut at selv om den bruker vilkårlige tall til å integrere så vil resultatene se veldig like ut for hver gang vi kjører koden. Vi fant ut at dette ikke er tilfellet med brute force Monte Carlo metoden, og at den vil gi forskjellige tall for hver gang vi kjører. Vi fant også ut at med litt arbeid så blir Gauss Laguerre en god metode, men Legendre er ikke noe spesielt bra for dette integralet.

Sist men ikke minst har vi laget et GitHub repository som inneholder kildekoden som vi har utviklet i dette projektet, i tillegg til et par håndplukkede resultater som skal brukte til å forsøke å vise at koden virker som den skal. GitHub repositoriet skal også vise en ganske detaljert logg av utviklingen av koden. GitHub repositoriet finner du her:

https://github.com/Erlendak/FYS3150_Project_3

Hvis du ønsker å clone reposetory kan du bruke denne lenken:

https://github.com/Erlendak/FYS3150_Project_3.git
