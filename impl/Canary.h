#ifndef CANARY
#define CANARY

#include <stdio.h>

#define CANARY_STR_IMPL_(x) #x
#define CANARY_STR(x) CANARY_STR_IMPL_(x)

#define CANARY_IS_INDEXABLE(arg) (sizeof(arg[0]))
#define CANARY_IS_ARRAY(arg) (IS_INDEXABLE(arg) && \
    (((void *) &arg) == ((void *) arg)))

#define CANARY_PRIMITIVE_DEFAULT                (0)
/* CHAR */
#define CANARY_PRIMITIVE_CHAR                   (1)
#define CANARY_PRIMITIVE_SIGNED_CHAR            (2)
#define CANARY_PRIMITIVE_UNSIGNED_CHAR          (3)
/* SHORT */
#define CANARY_PRIMITIVE_SHORT                  (4)
#define CANARY_PRIMITIVE_SHORT_INT              (5)
#define CANARY_PRIMITIVE_SIGNED_SHORT           (6)
#define CANARY_PRIMITIVE_SIGNED_SHORT_INT       (7)
#define CANARY_PRIMITIVE_UNSIGNED_SHORT         (8)
#define CANARY_PRIMITIVE_UNSIGNED_SHORT_INT     (9)
/* INT */
#define CANARY_PRIMITIVE_INT                    (10)
#define CANARY_PRIMITIVE_SIGNED                 (11)
#define CANARY_PRIMITIVE_SIGNED_INT             (12)
#define CANARY_PRIMITIVE_UNSIGNED               (13)
#define CANARY_PRIMITIVE_UNSIGNED_INT           (14)
/* LONG */
#define CANARY_PRIMITIVE_LONG                   (15)
#define CANARY_PRIMITIVE_LONG_INT               (16)
#define CANARY_PRIMITIVE_SIGNED_LONG            (17)
#define CANARY_PRIMITIVE_SIGNED_LONG_INT        (18)
/* LONG LONG */
#define CANARY_PRIMITIVE_LONG_LONG              (19)
#define CANARY_PRIMITIVE_LONG_LONG_INT          (20)
#define CANARY_PRIMITIVE_SIGNED_LONG_LONG       (21)
#define CANARY_PRIMITIVE_SIGNED_LONG_LONG_INT   (22)
#define CANARY_PRIMITIVE_UNSIGNED_LONG_LONG     (23)
#define CANARY_PRIMITIVE_UNSIGNED_LONG_LONG_INT (24)
/* FLOAT */
#define CANARY_PRIMITIVE_FLOAT                  (25)
/* DOUBLE */
#define CANARY_PRIMITIVE_DOUBLE                 (26)
/* LONG DOUBLE */
#define CANARY_PRIMITIVE_LONG_DOUBLE            (27)

#define CANARY_TYPECHECK(T)                         \
(_Generic((T),                                      \
    char:CANARY_PRIMITIVE_CHAR,                    \
    signed char:CANARY_PRIMITIVE_SIGNED_CHAR,      \
    unsigned char:CANARY_PRIMITIVE_UNSIGNED_CHAR,  \
    short:CANARY_PRIMITIVE_SHORT,                  \
    int:CANARY_PRIMITIVE_INT,                      \
    long:CANARY_PRIMITIVE_LONG,                    \
    float:CANARY_PRIMITIVE_FLOAT,                  \
    default:CANARY_PRIMITIVE_DEFAULT               \
))

#define CANARY_PRIMITIVE_CHAR_PLACEHOLDER                   "%c"
#define CANARY_PRIMITIVE_SIGNED_CHAR_PLACEHOLDER            "%c"
#define CANARY_PRIMITIVE_UNSIGNED_CHAR_PLACEHOLDER          "%c"
#define CANARY_PRIMITIVE_SHORT_PLACEHOLDER                  "%hi"
#define CANARY_PRIMITIVE_SHORT_INT_PLACEHOLDER              "%hi"
#define CANARY_PRIMITIVE_SIGNED_SHORT_PLACEHOLDER           "%hi"
#define CANARY_PRIMITIVE_SIGNED_SHORT_INT_PLACEHOLDER       "%hi"
#define CANARY_PRIMITIVE_UNSIGNED_SHORT_PLACEHOLDER         "%hu"
#define CANARY_PRIMITIVE_UNSIGNED_SHORT_INT_PLACEHOLDER     "%hu"
#define CANARY_PRIMITIVE_INT_PLACEHOLDER                    "%i"
#define CANARY_PRIMITIVE_SIGNED_PLACEHOLDER                 "%i"
#define CANARY_PRIMITIVE_SIGNED_INT_PLACEHOLDER             "%i"
#define CANARY_PRIMITIVE_UNSIGNED_PLACEHOLDER               "%u"
#define CANARY_PRIMITIVE_UNSIGNED_INT_PLACEHOLDER           "%u"
#define CANARY_PRIMITIVE_LONG_PLACEHOLDER                   "%li"
#define CANARY_PRIMITIVE_LONG_INT_PLACEHOLDER               "%li"
#define CANARY_PRIMITIVE_SIGNED_LONG_PLACEHOLDER            "%li"
#define CANARY_PRIMITIVE_SIGNED_LONG_INT_PLACEHOLDER        "%li"
#define CANARY_PRIMITIVE_LONG_LONG_PLACEHOLDER              "%lli"
#define CANARY_PRIMITIVE_LONG_LONG_INT_PLACEHOLDER          "%lli"
#define CANARY_PRIMITIVE_SIGNED_LONG_LONG_PLACEHOLDER       "%lli"
#define CANARY_PRIMITIVE_SIGNED_LONG_LONG_INT_PLACEHOLDER   "%lli"
#define CANARY_PRIMITIVE_UNSIGNED_LONG_LONG_PLACEHOLDER     "%llu"
#define CANARY_PRIMITIVE_UNSIGNED_LONG_LONG_INT_PLACEHOLDER "%llu"
#define CANARY_PRIMITIVE_FLOAT_PLACEHOLDER                  "%f"
#define CANARY_PRIMITIVE_DOUBLE_PLACEHOLDER                 "%lf"
#define CANARY_PRIMITIVE_LONG_DOUBLE_PLACEHOLDER            "%Lf"

char *PrimitivePlaceholder(int type) {
    switch (type) {
        case CANARY_PRIMITIVE_CHAR:
            return CANARY_PRIMITIVE_CHAR_PLACEHOLDER;
        case CANARY_PRIMITIVE_SIGNED_CHAR:
            return CANARY_PRIMITIVE_SIGNED_CHAR_PLACEHOLDER;
        case CANARY_PRIMITIVE_UNSIGNED_CHAR:
            return CANARY_PRIMITIVE_UNSIGNED_CHAR_PLACEHOLDER;
        case CANARY_PRIMITIVE_SHORT:
            return CANARY_PRIMITIVE_SHORT_PLACEHOLDER;
        case CANARY_PRIMITIVE_SHORT_INT:
            return CANARY_PRIMITIVE_SHORT_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_SIGNED_SHORT:
            return CANARY_PRIMITIVE_SIGNED_SHORT_PLACEHOLDER;
        case CANARY_PRIMITIVE_SIGNED_SHORT_INT:
            return CANARY_PRIMITIVE_SIGNED_SHORT_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_UNSIGNED_SHORT:
            return CANARY_PRIMITIVE_UNSIGNED_SHORT_PLACEHOLDER;
        case CANARY_PRIMITIVE_UNSIGNED_SHORT_INT:
            return CANARY_PRIMITIVE_UNSIGNED_SHORT_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_INT:
            return CANARY_PRIMITIVE_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_SIGNED:
            return CANARY_PRIMITIVE_SIGNED_PLACEHOLDER;
        case CANARY_PRIMITIVE_SIGNED_INT:
            return CANARY_PRIMITIVE_SIGNED_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_UNSIGNED:
            return CANARY_PRIMITIVE_UNSIGNED_PLACEHOLDER;
        case CANARY_PRIMITIVE_UNSIGNED_INT:
            return CANARY_PRIMITIVE_UNSIGNED_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_LONG:
            return CANARY_PRIMITIVE_LONG_PLACEHOLDER;
        case CANARY_PRIMITIVE_LONG_INT:
            return CANARY_PRIMITIVE_LONG_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_SIGNED_LONG:
            return CANARY_PRIMITIVE_SIGNED_LONG_PLACEHOLDER;
        case CANARY_PRIMITIVE_SIGNED_LONG_INT:
            return CANARY_PRIMITIVE_SIGNED_LONG_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_LONG_LONG:
            return CANARY_PRIMITIVE_LONG_LONG_PLACEHOLDER;
        case CANARY_PRIMITIVE_LONG_LONG_INT:
            return CANARY_PRIMITIVE_LONG_LONG_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_SIGNED_LONG_LONG:
            return CANARY_PRIMITIVE_SIGNED_LONG_LONG_PLACEHOLDER;
        case CANARY_PRIMITIVE_SIGNED_LONG_LONG_INT:
            return CANARY_PRIMITIVE_SIGNED_LONG_LONG_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_UNSIGNED_LONG_LONG:
            return CANARY_PRIMITIVE_UNSIGNED_LONG_LONG_PLACEHOLDER;
        case CANARY_PRIMITIVE_UNSIGNED_LONG_LONG_INT:
            return CANARY_PRIMITIVE_UNSIGNED_LONG_LONG_INT_PLACEHOLDER;
        case CANARY_PRIMITIVE_FLOAT:
            return CANARY_PRIMITIVE_FLOAT_PLACEHOLDER;
        case CANARY_PRIMITIVE_DOUBLE:
            return CANARY_PRIMITIVE_DOUBLE_PLACEHOLDER;
        case CANARY_PRIMITIVE_LONG_DOUBLE:
            return CANARY_PRIMITIVE_LONG_DOUBLE_PLACEHOLDER;
    }
    return CANARY_PRIMITIVE_DEFAULT;
}

#define CANARY_ACT(ACT) \
    { printf("BeginTest=%s\n", __func__); } \
    ACT; \
    do { printf("EndTest=%s\n", __func__); } while(0)

#define CANARY_TWEET_LOCATION(l)\
do { printf("Location="CANARY_STR(l)"\n"); } while(0)

#define CANARY_TWEET_BEGIN_UNIT(UNIT)\
do { printf("BeginUnit="CANARY_STR(UNIT)"\n"); } while(0)

#define CANARY_TWEET_END_UNIT(UNIT)\
do { printf("EndUnit="CANARY_STR(UNIT)"\n"); } while(0)

#define CANARY_TWEET_PRIMITIVE(VAR)                             \
do {                                                            \
    int type__canary_primitive__ = CANARY_TYPECHECK(VAR);       \
    printf("ID=%s\n", CANARY_STR(VAR));                         \
    printf("State=");                                           \
    printf(PrimitivePlaceholder(type__canary_primitive__), VAR);\
    printf("\n");                                               \
} while(0)

#define CANARY_TWEET_ERROR_STATE() \
do { printf("ERROR STATE\n"); } while (0)
#endif