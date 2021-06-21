# Create Project ----------------------------------------------------------


# Use R as calculator (type into console) ---------------------------------
2 + 2                  
7 - 8
4 * 6
144 / 12
8 ^ 9
2 + 2 * (2 ^ 6)

# Demonstrate the up key --------------------------------------------------

# Demonstrate the script --------------------------------------------------
8 * 6

# Creating objects --------------------------------------------------------
a <- 10 * 6
a

A <- 2 + 6 # CASE SENSISTIVE

A
a

# Overwriting objects -----------------------------------------------------
a <- 0
a

# Naming objects ----------------------------------------------------------
1object <- 3
!object <- 3
-object <- 3
object1 <- 3
object! <- 3
my.object <- 3
my_object <- 3
myObject <- 3

# Removing objects --------------------------------------------------------
rm(A)

# Data classes ------------------------------------------------------------
12.6         # NUMERIC
"cytometry"  # CHARACTER
TRUE         # LOGICAL

# Data structures ---------------------------------------------------------
# VECTOR
c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

c(1:10)

c("Mon", "Tue", "Wed", "Thur", "Fri")

c(Mon, Tue, Wed, Thur, Fri)           # CHARACTER STRINGS MUST BE IN QUOTES

c(FALSE, TRUE, TRUE, FALSE, FALSE)

c(false, true)                        # LOGICAL MUST BE IN CAPITALS

myVector <- c(F, F, F, T, T)          # BUT CAN BE ABBREVIATED

myVector

numbers <- c(1:10)

numbers * 2                           # OPERATIONS CAN BE PERFORMED ON VECTORS

numbers + 3

numbers + numbers


# LIST
list(1, 2, 3, "hello", FALSE)
list(myVector, 1, 2, 3, "hello", FALSE) # LISTS CAN CONTAIN LISTS AND VECTORS
list(myVector, c(1, 2, 3), "hello", FALSE)

# DATA.FRAME (COME BACK TO)

# Subsetting vectors ------------------------------------------------------
days <- c("Mon", "Tue", "Wed", "Thur", "Fri")
days
days[1]
days[4]
days[c(1, 3, 4)]
days[1:4]
days[-5]

# Subsetting lists --------------------------------------------------------
my_list <- list(a = 1:3, b = TRUE, c = "laser")
my_list[[1]]
my_list[["a"]]
my_list[["a"]][2]

# Functions ---------------------------------------------------------------
myValues <- c(1:100)
myValues
mean(myValues)
median(myValues)
min(myValues)
max(myValues)
sum(myValues)
sd(myValues)
class(myValues)
length(myValues) # SOME FUNCTIONS RETURN SINGLE VALUES (AGGREGATE FUNCTIONS)
log(myValues)    # OTHERS RETURN A VALUE FOR EACH COMPONENT OF THE VECTOR
log10(myValues)  # CAREFUL: DIFFERENCE BETWEEN LOG10 AND LOG
mySqrt <- sqrt(myValues)
mySqrt
?rnorm           # HELP ON HOW TO USE A FUNCTION
rnorm(100, mean = 5)
hist(rnorm(100, mean = 5))

# Data frames -------------------------------------------------------------
id <- 1:200
group <- c(rep("Vehicle", 100), rep("Drug", 100))
response <- c(rnorm(100, mean = 25, sd = 5), rnorm(100, mean = 23, sd = 5))
myData <- data.frame(Patient = id, Treatment = group, Response = response)
myData # CTRL+L
head(myData) # REMIND DATA WILL BE DIFFERENT BECAUSE OF RANDOM SEED
head(myData, 12)
tail(myData, 10)
dim(myData)
str(myData)
summary(myData)

# SUBSETTING DATA.FRAMES ----
myData[1, 2] #[ROWS, COLUMNS]
myData[2, 3]
myData[1:20, 2:3]
myData[1:20, ]
myData[, 3]
myData[, "Response"]
myData$Response

responders <- myData$Response > 26
responders
myData[responders, ]
myData[myData$Treatment == "Vehicle" & myData$Response <= 23, ]
myData[myData$Treatment == "Vehicle" | myData$Response >= 21, ]
myData[myData$Treatment != "Vehicle" | myData$Response > 24, ]
age <- round(rnorm(200, mean = 40, sd = 20))
myData$Age <- age
head(myData)
 

# Reading data into R -----------------------------------------------------
pokemon <- read.csv("data/Pokemon.csv", header = T) # EXPLAIN HEADER (DEFAULT IS TRUE)
dim(pokemon)
head(pokemon)
str(pokemon)
summary(pokemon)

# Installing packages -----------------------------------------------------
install.packages("dplyr")

# Loading packages --------------------------------------------------------
library(dplyr)

# Selecting columns -------------------------------------------------------
select(pokemon, Pokemon, Def, Atk)
select(pokemon, 2, 5, 4)

# Filter rows -------------------------------------------------------------
filter(pokemon, HP > 4)
filter(pokemon, HP > 4, Type.I == "Bug")

# Summarising a group -----------------------------------------------------
grouped <- group_by(pokemon, Type.I)

summarise(grouped, "Median_Atk" = median(Atk))

# Mutating new columns ----------------------------------------------------
mutate(pokemon, Atk_per_Def = Atk / Def)
mutate(pokemon, 
       Atk_per_Def = Atk / Def,
       HP_min_Atk = HP - Atk)

# Arranging rows ----------------------------------------------------------
arrange(pokemon, Total)
arrange(pokemon, -Total)

# Chaining operations together with pipe ----------------------------------
pokemon %>% summary()

pokemon %>%
    filter(Type.I != "Grass") %>%
    mutate(Atk_per_Def = Atk / Def) %>%
    group_by(Type.II) %>% 
    summarise(Atk_per_Def = median(Atk_per_Def)) %>%
    arrange(Atk_per_Def)

