# Regenerate the Boucot palaeolithology dataset from the source CSV tables

The `boucot_paleolithology_combined` shapefile in `gprm/Data` had five text columns
(`Lithology`, `Formation`, `LithComm`, `Period`, `AgeComm`) destroyed by a numeric type
coercion during an earlier CSV-to-shapefile conversion — every value is `0.0` or null,
including the free-text lithology description that is the dataset's most informative field.
Rather than ship that, we rebuild the bundled dataset by joining the 28 pristine source
CSVs from Boucot, Chen & Scotese (2013) back onto the shapefile's GPlates attributes
(`PLATEID1`, `FROMAGE`, `TOAGE`, `FEATURE_ID`), which the CSVs do not carry.

## Considered Options

Shipping the shapefile as-is was rejected because `Lithology.value_counts()` returning
`{0.0: 7586}` reads as data rather than corruption. Dropping the five columns was rejected
because it makes the loss permanent and silent.

## Consequences

The join key is `(time slice, LithCode, LithID)` — `LithID` restarts within each of the 28
map files, so it is not unique on its own. The time slice is recovered from `FROMAGE`,
whose 28 distinct values map 1:1 onto the source files. Validated at 8698/8698 rows matched,
with coordinates agreeing to better than 0.001° for all but 15 rows and 5 duplicate keys
needing hand-reconciliation.

The source CSVs live outside the repository and are cp1252-encoded with inconsistent
headers across files (`LithologyIDNumber` vs `LithNumber`, singular vs plural `*Comment(s)`,
a leading space on `LithologyCode` in `Map07 Eifel v5.csv`, stray unnamed columns). The
build script normalises these; the build is therefore not reproducible from a clone alone,
which is the accepted cost of not redistributing SEPM's tables.
