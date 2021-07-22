#include "jstools.hh"
#include "out.hh"

void jstools::assertKey(const json11::Json & json, const std::string & key)
{
    if (json.object_items().count(key) == 0)
    {
        out("Config json does not contain required key:", key);
        exit(101);
    }
}

bool jstools::contains(const json11::Json & json, const std::string & key)
{
    if (!json.is_object()) return false;
    return (json.object_items().count(key) > 0);
}

void jstools::readInt(const json11::Json & json, const std::string & key, int & var)
{
    assertKey(json, key);
    if (!json[key].is_number())
    {
        out(key, "is not a number");
        exit(102);
    }
    var = json[key].int_value();
    out(key, var);
}

void jstools::readDouble(const json11::Json & json, const std::string & key, double & var)
{
    assertKey(json, key);
    if (!json[key].is_number())
    {
        out(key, "is not a number");
        exit(102);
    }
    var = json[key].number_value();
    out(key, var);
}

void jstools::readBool(const json11::Json & json, const std::string & key, bool & var)
{
    assertKey(json, key);
    if (!json[key].is_bool())
    {
        out(key, "is not a bool");
        exit(102);
    }
    var = json[key].bool_value();
    out(key, (var ? "true" : "false"));
}

void jstools::readArray(const json11::Json &json, const std::string &key, json11::Json::array & var)
{
    assertKey(json, key);
    if (!json[key].is_array())
    {
        out(key, "is not an array");
        exit(102);
    }
    var = json[key].array_items();
}

void jstools::readString(const json11::Json &json, const std::string &key, std::string & var)
{
    assertKey(json, key);
    if (!json[key].is_string())
    {
        out(key, "is not a string");
        exit(102);
    }
    var = json[key].string_value();
    out(key, var);
}

void jstools::readObject(const json11::Json &json, const std::string &key, json11::Json::object &var)
{
    assertKey(json, key);
    if (!json[key].is_object())
    {
        out(key, "is not an object");
        exit(102);
    }
    var = json[key].object_items();
}
